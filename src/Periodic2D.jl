@with_kw struct Periodic2D_params
  # Number of MPI processes
  nprocs::Int = 1
  rank::Int = 1

  # Case name
  case::String = "test"
  iter::Int = 0

  # Geometric parameters
  H::Float64 = 1.1 # Height of the domain
  Lf::Float64 = 1.0 # Length of the entire domain
  hs::Float64 = 0.01 # Floating structure thickness

  # Wave Parameters
  kλ::Float64 = 3.0 # Wave number
  η₀::Float64 = 0.01 # surface elevation

  # Time Numerics
  ρ∞::Float64 = 0.5
  dt::Float64 = 0.1
  t0::Float64 = 0.0
  tF::Float64 = 0.1

  # Physical Parameters
  ρf::Float64 = 1000.0 # Fluid density
  ρs::Float64 = 100.0 # Solid density
  g::Float64 = 9.81 # Acceleration due to gravity

  # Robin parameter
  αf::Float64 = (
    ρs * hs * (1 - ((2 * ρ∞ - 1) / (1 + ρ∞)))
  ) / (
    dt *
    (1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞)) *
    (1 - (ρ∞ / (1 + ρ∞)))
  )
  αs::Float64 = 0.0

  # Postprocessing Parameters
  vtkoutput::Bool = false
end

function periodic2D(distribute, parts, params::Periodic2D_params)
  # Initializing the timer
  to = TimerOutput()

  @timeit to "Model & Setup" begin
    # ranks
    ranks = distribute(LinearIndices((prod(parts),)))

    # Unpack the parameters
    @unpack_Periodic2D_params params

    # Wave parameters
    ω = sqrt(g * kλ * tanh(kλ * H)) # Wave frequency in radians
    ϕ = 0 # wave phase difference

    # Derived parameters
    T = ((ω * ω) / (kλ * kλ)) * (
      (ρf) / (kλ * tanh(kλ * H)) + ρs * hs
    )

    # Creating the mesh

    # Defining the model
    domain = (0, Lf, 0, H)
    partition = (trunc(Lf * 7), trunc(H * 7))

    function f_z(x)
      if x == H
        return H
      end
      i = x / (H / trunc(H * 7))
      return H - H / (1.2^i)
    end
    map_cord(x) = VectorValue(x[1], f_z(x[2]))

    model = UnstructuredDiscreteModel(
      CartesianDiscreteModel(
        ranks,
        parts,
        domain,
        partition,
        map = map_cord,
        isperiodic = (true, false),
      ),
    )

    # Define tags in the model
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "FloatingSolid", [6])
    add_tag_from_tags!(labels, "Bed", [1, 5, 2])
    add_tag_from_tags!(labels, "LeftPoint", [3])
    add_tag_from_tags!(labels, "RightPoint", [4])
    Geometry.add_tag_from_tags_complementary!(
      labels,
      "!!FloatingSolid",
      ["FloatingSolid"],
    )
    Geometry.add_tag_from_tags_setdiff!(
      labels,
      "!FloatingSolid",
      ["!!FloatingSolid"],
      ["LeftPoint", "RightPoint"],
    )

    # Define reference FE (P2/P1 pair)
    fe_order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, fe_order)
    reffeₚ = ReferenceFE(lagrangian, Float64, fe_order - 1)
    reffeₛ = ReferenceFE(lagrangian, Float64, fe_order)

    # Define triangulation and integration measure
    degree = 2 * fe_order + 1

    Ωs = get_triangulations(model, "FloatingSolid")
    dΩs = Measure(Ωs, degree)

    Ωf = Triangulation(model)
    dΩf = Measure(Ωf, degree)

    Σs = BoundaryTriangulation(model, tags = ["FloatingSolid"])
    dΣs = Measure(Σs, degree)

    # Define Test FE Functions
    mfs = BlockMultiFieldStyle(2, (1, 2), (1, 2, 3))
    V = TestFESpace(
      Ωf,
      reffeᵤ,
      dirichlet_tags = ["Bed"],
      dirichlet_masks = [(false, true)],
    )
    Q = TestFESpace(Ωf, reffeₚ, conformity = :H1)
    S = TestFESpace(Ωs, reffeₛ, dirichlet_tags = ["!FloatingSolid"])
    Y = MultiFieldFESpace([S, V, Q]; style = mfs)

    # Define pressure boundary condition functions
    ux(t, x) = (
      ω *
      η₀ *
      (cosh(kλ * x[2]) / sinh(kλ * H)) *
      cos(kλ * x[1] - ω * t + ϕ)
    )
    uy(t, x) = (
      ω *
      η₀ *
      (sinh(kλ * x[2]) / sinh(kλ * H)) *
      sin(kλ * x[1] - ω * t + ϕ)
    )
    u_field(t) = x -> VectorValue(ux(t, x), uy(t, x))
    pres_field(t) = x -> (
      ρf * g * (H - x[2]) +
      ρf *
      ((η₀ * ω * ω) / kλ) *
      (cosh(kλ * x[2]) / sinh(kλ * H)) *
      cos(kλ * x[1] - ω * t + ϕ)
    )
    float_solid(t) = x -> η₀ * cos(kλ * x[1] - ω * t + ϕ)
    zero_scalar(t) = x -> 0.0

    # Dual use (vel/acc) terms
    pressure_vel_acc(t) = x -> 0.0 # Initial vel/acc pressure
    velocity_vel_acc(t) = x -> VectorValue(0.0, 0.0)
    sol_dis_vel_acc(t) = x -> 0.0 # Initial vel/acc solid displacement

    # Define Trial FE Functions
    U = TransientTrialFESpace(V, [u_field]) # Trial fluid velocity
    P = TransientTrialFESpace(Q) # Trial Pressure
    D = TransientTrialFESpace(S, [zero_scalar]) # Trial membrane displacement
    X = MultiFieldFESpace([D, U, P]; style = mfs) # Trial Multifield

    # Damping function
    α(x) = 1.0
    αvec(x) = VectorValue(α(x), α(x))
    alpha = interpolate_everywhere([α, αvec, α], X(0.0))
    x_base(t) = interpolate_everywhere(
      [float_solid(t), u_field(t), pres_field(t)],
      X(t),
    )

    nᵥ = VectorValue(0.0, 1.0) # Along the Y direction
    nᵤ = VectorValue(1.0, 0.0) # Along the X direction

    # Weak form definition
    jac(t, (dd, du, dp), (s, v, q)) = (
      ∫(-(∇ ⋅ v) * dp)dΩf +
      ∫((∇ ⋅ du) * q)dΩf +
      ∫(T * ((∇(dd) ⋅ nᵤ) * (∇(v ⋅ nᵥ) ⋅ nᵤ)))dΣs +
      ∫(ρf * g * dd * (v ⋅ nᵥ))dΣs +
      ∫(αf * (du ⋅ nᵥ) * (v ⋅ nᵥ))dΣs +
      ∫(T * (∇(dd) ⋅ nᵤ) * (∇(s) ⋅ nᵤ))dΣs +
      ∫(ρf * g * dd * s)dΣs -
      ∫(dp * (∇(s) ⋅ nᵥ))dΩs
    )
    jac_t(t, (dtd, dtu, dtp), (s, v, q)) = (
      ∫(ρf * (v ⋅ dtu))dΩf -
      ∫(αf * dtd * (v ⋅ nᵥ))dΣs +
      ∫(ρf * (dtu ⋅ nᵥ) * s)dΩs
    )
    jac_tt(t, (dttd, dttu, dttp), (s, v, q)) = (
      ∫(ρs * hs * (v ⋅ nᵥ) * dttd)dΣs +
      ∫(ρs * hs * s * dttd)dΣs
    )
    l(t, (s, v, q)) = (
      ∫((v ⋅ nᵥ) * -ρf * g + q * 0.0)dΩf +
      ∫(s * -ρf * g)dΩs
    )

    # Build affine FE operator
    assem = SparseMatrixAssembler(
      SparseMatrixCSR{0,Float64,Int},
      Vector{Float64},
      X,
      Y,
    )
    op = TransientLinearFEOperator(
      (jac, jac_t, jac_tt),
      l,
      X,
      Y;
      assembler = assem,
    )

    # System Solver Definition
    sol_param_fluid_dir = datadir("periodic2d", "solver_parameters_fluid.xml")
    fluid_block = TrilinosSolve(sol_param_fluid_dir)

    sol_param_solid_dir = datadir("periodic2d", "solver_parameters_solid.xml")
    solid_block = TrilinosSolve(sol_param_solid_dir)

    coeffs = [
      1.0 1.0
      0.0 1.0
    ]

    bblocks = [
      LinearSystemBlock() LinearSystemBlock()
      LinearSystemBlock() LinearSystemBlock()
    ]

    prec = BlockTriangularSolver(
      bblocks,
      [solid_block, fluid_block],
      coeffs,
      :upper,
    )
    solver = FGMRESSolver(
      20,
      prec;
      Pl = nothing,
      restart = false,
      m_add = 1,
      maxiter = 1000,
      atol = 1e-8,
      rtol = 1.0e-7,
      verbose = i_am_main(ranks),
    )
    sys_solver = DiscreteDampingSolver(solver, alpha, x_base)

    # ODE Solver
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(sys_solver, dt, ρ∞)

    # Imposing initial conditions
    x_initial = interpolate_everywhere(
      [float_solid(t0), u_field(t0), pres_field(t0)],
      X(t0),
    )
    v_initial = interpolate_everywhere(
      [sol_dis_vel_acc(t0), velocity_vel_acc(t0), pressure_vel_acc(t0)],
      X(t0),
    )
    a_initial = interpolate_everywhere(
      [sol_dis_vel_acc(t0), velocity_vel_acc(t0), pressure_vel_acc(t0)],
      X(t0),
    )

    # Initializing the solving process
    xₜ = solve(solver_ode, op, t0, tF, (x_initial, v_initial, a_initial))

    @timeit to "Postprocessing Initial Values" begin
      @vtk postprocess_dir = mkpath(
        datadir("periodic2d", case, "postprocess"),
      )
      @vtk vtk_dir = mkpath("$(postprocess_dir)/tmp")

      # Creating pvd files
      @vtk pvd_Ωf = createpvd(
        Ωf,
        ranks,
        "$(postprocess_dir)/fluids_in_twoD",
      )
      @vtk pvd_Σs = createpvd(
        Σs,
        ranks,
        "$(postprocess_dir)/solids_in_twoD",
      )

      # Storing initial conditions in the pvd file
      @vtk pvd_Ωf[0] = createvtk(
        Ωf,
        "$(vtk_dir)/results_fluids_in_twoD_0",
        cellfields = ["u" => x_initial[2], "pressure" => x_initial[3]],
      )
      @vtk pvd_Σs[0] = createvtk(
        Σs,
        "$(vtk_dir)/results_solids_in_twoD_0",
        cellfields = ["displacement" => x_initial[1]],
      )
    end # Ending the postprocessing initial values timer
  end # Ending the model and setup timer

  timestep = 0
  n_timesteps = round(Int, (tF - t0) / dt)
  @timeit to "Time Stepping Loop Aggregate" begin
    solve_timer = begin_timed_section!(to, "First solve")
    # Solving and storing results in the pvd file
    with_logger(SimpleLogger(stderr, Logging.Error)) do
      for (tn, (dhn, uhn, phn)) in xₜ
        @vtk pvd_Ωf[tn] = createvtk(
          Ωf,
          "$(vtk_dir)/results_fluids_in_twoD_$tn",
          cellfields = ["u" => uhn, "pressure" => phn],
        )
        @vtk pvd_Σs[tn] = createvtk(
          Σs,
          "$(vtk_dir)/results_solids_in_twoD_$tn",
          cellfields = ["displacement" => dhn],
        )
        if solve_timer !== nothing
          end_timed_section!(to, solve_timer)
          solve_timer = nothing
        end
        i_am_main(ranks) && println("$(timestep + 1) timestep solved")
        timestep += 1
        if timestep == n_timesteps - 1
          solve_timer = begin_timed_section!(to, "Subsequent Solves")
        end
      end
    end
    if solve_timer !== nothing
      end_timed_section!(to, solve_timer)
    end
  end # Ending the time stepping loop timer

  @timeit to "Saving pvds" begin
    # Saving the pvd files
    @vtk savepvd(pvd_Ωf)
    @vtk savepvd(pvd_Σs)
  end # Ending the saving results timer

  # Display the timer summary
  if i_am_main(ranks)
    println("\n--- Simulation Performance Summary ---")
    show(to)
    println("\n--------------------------------------")
  end

  # Return convergence and timer data for the solver at time tF
  outer_iter = solver.log.num_iters
  outer_iter_array = collect(0:outer_iter)
  outer_residuals = solver.log.residuals
  outer_residuals = outer_residuals[outer_iter_array .+ 1]

  return (
    fluid = (
      num_iters = fluid_block.log.num_iters,
      residual = fluid_block.log.residual,
      solve_time = fluid_block.log.solve_time,
    ),
    solid = (
      num_iters = solid_block.log.num_iters,
      residual = solid_block.log.residual,
      solve_time = solid_block.log.solve_time,
    ),
    outer = (
      iter_array = outer_iter_array,
      residuals = outer_residuals,
      timer = to,
    ),
  )
end
# Ending FSI Problem Definition
###############################

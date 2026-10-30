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
  ) + 2 * ρf * g * dt / ((1 + ρ∞) * (3 - ρ∞))
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
    # Same mixed spaces as WSI; the entire top boundary is a membrane.
    fe_order = 2
    reffeᵤ = ReferenceFE(raviart_thomas, Float64, fe_order)
    reffeₚ = ReferenceFE(lagrangian, Float64, fe_order)

    degree = 2 * fe_order + 2
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf, degree)

    Σs = BoundaryTriangulation(model, tags = ["FloatingSolid"])
    dΣs = Measure(Σs, degree)
    Λ = SkeletonTriangulation(Σs)
    dΛ = Measure(Λ, degree)
    nΛ = get_normal_vector(Λ)

    # The fluid block contains the multiplier, velocity, and pressure.
    mfs = BlockMultiFieldStyle(2, (1, 3), (1, 2, 3, 4))
    V = TestFESpace(model, reffeᵤ, conformity = :HDiv, dirichlet_tags = ["Bed"])
    Q = TestFESpace(Ωf, reffeₚ, conformity = :L2)
    E = TestFESpace(Σs, reffeₚ, conformity = :H1)
    L = TestFESpace(Σs, reffeₚ, conformity = :L2)
    Y = MultiFieldFESpace([E, L, V, Q]; style = mfs)

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

    zero_vector(t) = x -> VectorValue(0.0, 0.0)
    # The multiplier is the pressure in excess of hydrostatic restoring pressure.
    membrane_pressure(t) = x -> pres_field(t)(x) - ρf * g * float_solid(t)(x)

    # Periodicity is inherited from the model; solid and multiplier live on its surface.
    U = TransientTrialFESpace(V, [zero_vector])
    P = TransientTrialFESpace(Q)
    Et = TransientTrialFESpace(E)
    Lt = TransientTrialFESpace(L)
    X = MultiFieldFESpace([Et, Lt, U, P]; style = mfs)

    nᵥ = VectorValue(0.0, 1.0)

    # Dimensional WSI formulation, with no uncovered free surface or damping zones.
    jac(t, (dη, dλ, du, dp), (s, μ, v, q)) = (
      # Membrane dynamics
      ∫(T * (∇(dη) ⋅ ∇(s)))dΣs -
      ∫(dλ * s)dΣs +
      # Robin interface: kinematics and membrane traction
      ∫(αf * μ * (du ⋅ nᵥ))dΣs +
      ∫(T * (∇(dη) ⋅ ∇(μ)))dΣs -
      ∫(T * (jump(μ) * nΛ.⁺) ⋅ mean(∇(dη)))dΛ -
      ∫(dλ * μ)dΣs +
      # Fluid momentum and incompressibility
      ∫(-(∇ ⋅ v) * dp)dΩf +
      ∫((ρf * g * dη + dλ) * (v ⋅ nᵥ))dΣs +
      ∫((∇ ⋅ du) * q)dΩf
    )
    jac_t(t, (dtη, dtλ, dtu, dtp), (s, μ, v, q)) = (
      ∫(-αf * dtη * μ)dΣs +
      ∫(ρf * (v ⋅ dtu))dΩf
    )
    jac_tt(t, (dttη, dttλ, dttu, dttp), (s, μ, v, q)) = (
      ∫(ρs * hs * dttη * s)dΣs +
      ∫(ρs * hs * dttη * μ)dΣs
    )
    l(t, (s, μ, v, q)) = (
      ∫(-ρf * g * (v ⋅ nᵥ))dΩf +
      ∫(0.0 * s)dΣs + ∫(0.0 * μ)dΣs +
      ∫(0.0 * q)dΩf + ∫(0.0 * jump(μ))dΛ
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
    sys_solver = SolverWrapper(solver)

    # ODE Solver
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(sys_solver, dt, ρ∞)

    # Imposing initial conditions
    x_initial = interpolate_everywhere(
      [float_solid(t0), membrane_pressure(t0), u_field(t0), pres_field(t0)],
      X(t0),
    )
    v_initial = interpolate_everywhere(
      [zero_scalar(t0), zero_scalar(t0), zero_vector(t0), zero_scalar(t0)],
      X(t0),
    )
    a_initial = interpolate_everywhere(
      [zero_scalar(t0), zero_scalar(t0), zero_vector(t0), zero_scalar(t0)],
      X(t0),
    )

    # Initializing the solving process
    xₜ = solve(solver_ode, op, t0, tF, (x_initial, v_initial, a_initial))

    @timeit to "Postprocessing Initial Values" begin
      @vtk postprocess_dir = mkpath(
        datadir("periodic2d", case, "postprocess"),
      )
      @vtk vtk_dir = mkpath(joinpath(postprocess_dir, "tmp"))

      # Creating pvd files
      @vtk pvd_Ωf = createpvd(
        Ωf,
        ranks,
        joinpath(postprocess_dir, "fluids_in_twoD"),
      )
      @vtk pvd_Σs = createpvd(
        Σs,
        ranks,
        joinpath(postprocess_dir, "solids_in_twoD"),
      )

      # Storing initial conditions in the pvd file
      @vtk pvd_Ωf[t0] = createvtk(
        Ωf,
        joinpath(vtk_dir, "results_fluids_in_twoD_0"),
        cellfields = ["u" => x_initial[3], "pressure" => x_initial[4]],
      )
      @vtk pvd_Σs[t0] = createvtk(
        Σs,
        joinpath(vtk_dir, "results_solids_in_twoD_0"),
        cellfields = ["displacement" => x_initial[1], "membrane_pressure" => x_initial[2]],
      )
    end # Ending the postprocessing initial values timer
  end # Ending the model and setup timer

  timestep = 0
  @timeit to "Time Stepping Loop Aggregate" begin
    solve_timer = begin_timed_section!(to, "First solve")
    # Solving and storing results in the pvd file
    with_logger(SimpleLogger(stderr, Logging.Error)) do
      for (tn, (ηhn, λhn, uhn, phn)) in xₜ
        @vtk pvd_Ωf[tn] = createvtk(
          Ωf,
          joinpath(vtk_dir, "results_fluids_in_twoD_$tn"),
          cellfields = ["u" => uhn, "pressure" => phn],
        )
        @vtk pvd_Σs[tn] = createvtk(
          Σs,
          joinpath(vtk_dir, "results_solids_in_twoD_$tn"),
          cellfields = ["displacement" => ηhn, "membrane_pressure" => λhn],
        )
        end_timed_section!(to, solve_timer)
        i_am_main(ranks) && println("$(timestep + 1) timestep solved")
        timestep += 1
        solve_timer = begin_timed_section!(to, "Subsequent Solves")
      end
    end
    end_timed_section!(to, solve_timer)
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

  solver_stats = (
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
  GC.gc(true)
  return solver_stats
end
# Ending FSI Problem Definition
###############################

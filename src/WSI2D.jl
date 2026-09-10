@with_kw struct WSI2D_params
  # Number of MPI processes
  nprocs::Int = 1
  rank::Int = 1

  # Case name
  case::String = "test"

  # Geometrical Parameters
  H::Float64
  Lm::Float64
  Lf::Float64
  hs::Float64
  meshpath::String

  # Damping Parameters
  Lfd::Float64
  Lfd1::Float64
  Ld::Float64
  Ld1::Float64

  # Wave Parameters
  kλ::Float64
  η₀::Float64
  ϕ::Float64
  ω::Float64

  # Time Numerics
  ρ∞::Float64
  dt::Float64
  t0::Float64
  tF::Float64

  # Physical Parameters
  ρf::Float64
  ρs::Float64
  g::Float64
  T::Float64

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

function ramp(time; t_ramp = 2.0)
  return clamp(time / t_ramp, 0.0, 1.0)
end

#################################
# Starting FSI Problem Definition
function wsi2d(distribute, parts, params::WSI2D_params)

  # Initializing the timer
  to = TimerOutput()

  @timeit to "Model & Setup" begin

    # ranks
    ranks = distribute(LinearIndices((prod(parts),)))

    # Unpack the parameters
    @unpack_WSI2D_params params

    # Defining the model
    model = UnstructuredDiscreteModel(
      GmshDiscreteModel(ranks, meshpath, renumber = false),
    )
    labels = get_face_labeling(model)
    Geometry.add_tag_from_tags_complementary!(
      labels,
      "!!FreeSurface",
      ["FreeSurface"],
    )
    Geometry.add_tag_from_tags_setdiff!(
      labels,
      "!FreeSurface",
      ["!!FreeSurface"],
      ["LeftPoint", "RightPoint"],
    )
    Geometry.add_tag_from_tags_complementary!(
      labels,
      "!FloatingSolid",
      ["FloatingSolid"],
    )

    # Define reference FE (P2/P1 taylor hood elements for fe_order = 2)
    fe_order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, fe_order)
    reffeₚ = ReferenceFE(lagrangian, Float64, fe_order - 1)
    reffeₛ = ReferenceFE(lagrangian, Float64, fe_order)

    # Define triangulation and integration measure
    degree = 2 * fe_order + 1
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf, degree)

    Ωfs, Ωs = get_triangulations(model, "FreeSurface", "FloatingSolid")

    dΩfs = Measure(Ωfs, degree)

    dΩs = Measure(Ωs, degree)

    Σs = BoundaryTriangulation(model, tags = ["FloatingSolid"])
    dΣs = Measure(Σs, degree)

    Σfs = BoundaryTriangulation(model, tags = ["FreeSurface"])
    dΣfs = Measure(Σfs, degree)

    Σinlet = BoundaryTriangulation(model, tags = ["Inlet"])
    dΣinlet = Measure(Σinlet, degree)

    Σoutlet = BoundaryTriangulation(model, tags = ["Outlet"])
    dΣoutlet = Measure(Σoutlet, degree)

    # Define Test FE Functions
    mfs = BlockMultiFieldStyle(3, (1, 1, 2), (1, 2, 3, 4))
    V = TestFESpace(
      Ωf,
      reffeᵤ,
      dirichlet_tags = ["Bed"],
      dirichlet_masks = [(false, true)],
    ) # Test function for fluid velocity
    Q = TestFESpace(Ωf, reffeₚ) # Test function for Pressure
    S = TestFESpace(
      Ωs,
      reffeₛ,
      dirichlet_tags = ["!FloatingSolid"],
    ) # Test function for Solid Displacement
    Sfs = TestFESpace(
      Ωfs,
      reffeₛ,
      dirichlet_tags = ["!FreeSurface"],
    ) # Test function for free surface elevation
    Y = MultiFieldFESpace([S, Sfs, V, Q]; style = mfs) # Test Multifield

    Gs = select_triangulation(Ωs, Ωfs)
    dGs = Measure(Gs, degree)
    Gfs = select_triangulation(Ωfs, Ωs)
    dGfs = Measure(Gfs, degree)
    ns = get_normal_vector(Gs)
    nfs = get_normal_vector(Gfs)

    # Defining variable fields
    # Displacement terms
    function u_field(t)
      ux(t, x) = (
        -ω *
        η₀ *
        (cosh(kλ * x[2]) / sinh(kλ * H)) *
        cos(kλ * x[1] - ω * t + ϕ)
      )
      uy(t, x) = (
        -ω *
        η₀ *
        (sinh(kλ * x[2]) / sinh(kλ * H)) *
        sin(kλ * x[1] - ω * t + ϕ)
      )
      function inner_function(x)
        if x[1] <= Lfd
          ans = VectorValue(ux(t, x), uy(t, x)) * ramp(t)
        else
          ans = VectorValue(0.0, 0.0)
        end
        return ans
      end
      return inner_function
    end

    function pres_field(t)
      function inner_function(x)
        if x[1] <= Lfd
          ans = (
            ρf * g * (H - x[2]) -
            ((ρf * η₀ * ω * ω) / kλ) *
            (cosh(kλ * x[2]) / sinh(kλ * H)) *
            cos(kλ * x[1] - ω * t + ϕ) *
            ramp(t)
          )
        else
          ans = ρf * g * (H - x[2])
        end
        return ans
      end
      return inner_function
    end

    function free_surface_field(t)
      function inner_function(x)
        if x[1] <= Lfd
          ans = (
            -1.0 *
            ((η₀ * ω * ω) / (kλ * g)) *
            (cosh(kλ * x[2]) / sinh(kλ * H)) *
            cos(kλ * x[1] - ω * t + ϕ) *
            ramp(t)
          )
        else
          ans = 0.0
        end
        return ans
      end
      return inner_function
    end

    float_solid(t) = x -> 0.0 # Explicit condition for solid displacement

    # Dual use (vel/acc) terms
    pressure_vel_acc(t) = x -> 0.0 # Initial vel/acc terms for pressure
    velocity_vel_acc(t) = x -> VectorValue(0.0, 0.0) # Initial vel/acc velocity
    sol_dis_vel_acc(t) = x -> 0.0 # Initial vel/acc solid displacement
    fs_ele_vel_acc(t) = x -> 0.0 # Initial vel/acc free surface elevation

    # Define Trial FE Functions
    U = TransientTrialFESpace(V, [u_field]) # Trial function for fluid velocity
    P = TransientTrialFESpace(Q) # Trial function for Pressure
    D = TransientTrialFESpace(S, [float_solid]) # Trial Solid Displacement
    Dfs = TransientTrialFESpace(Sfs, [float_solid]) # Trial free surface
    X = MultiFieldFESpace([D, Dfs, U, P]; style = mfs) # Trial Multifield

    # Damping function
    function α(x)
      c1 = 0.1
      c2 = 10.0
      function frac(x, Ld, Ld1, Lfd, Lfd1)
        if x[1] ≈ Ld1 || x[1] > Ld1
          return 1.0 # 100% outlet damping after Ld1
        elseif x[1] >= Ld
          return (x[1] - Ld) / (Ld1 - Ld)
        elseif x[1] ≈ Lfd1 || x[1] < Lfd1
          return 1.0 # 100% inlet damping before Lfd1
        elseif x[1] <= Lfd
          return (Lfd - x[1]) / (Lfd - Lfd1)
        else
          return 0.0 # No damping in the middle region
        end
      end
      frac_x = frac(x, Ld, Ld1, Lfd, Lfd1)
      return (
        (1.0 - c1 * frac_x * frac_x) *
        (1.0 - (1.0 - exp(c2 * frac_x * frac_x)) / (1.0 - exp(c2)))
      )
    end

    αvec(x) = VectorValue(α(x), α(x))
    alpha = interpolate_everywhere([α, α, αvec, α], X(0.0))
    x_base(t) = interpolate_everywhere(
      [float_solid(t), free_surface_field(t), u_field(t), pres_field(t)],
      X(t),
    )

    nᵤ = VectorValue(1.0, 0.0)
    nᵥ = VectorValue(0.0, 1.0)

    # Weak form definition
    jac(t, (dd, dη, du, dp), (s, γ, v, q)) = (
      ∫(-(∇ ⋅ v) * dp)dΩf +
      ∫((∇ ⋅ du) * q)dΩf +
      ∫(T * ((∇(dd) ⋅ nᵤ) * (∇(v ⋅ nᵥ) ⋅ nᵤ)))dΣs +
      ∫(ρf * g * dd * (v ⋅ nᵥ))dΣs +
      ∫(ρf * g * dη * (v ⋅ nᵥ))dΣfs +
      ∫(αf * (du ⋅ nᵥ) * (v ⋅ nᵥ))dΣs +
      ∫(αf * (du ⋅ nᵥ) * (v ⋅ nᵥ))dΣfs +
      ∫(T * (∇(dd) ⋅ nᵤ) * (∇(s) ⋅ nᵤ))dΣs +
      ∫(ρf * g * dd * s)dΣs -
      ∫(dp * (∇(s) ⋅ nᵥ))dΩs +
      ∫(dp * s * (ns ⋅ nᵥ))dGs +
      ∫(ρf * g * dη * γ)dΣfs -
      ∫(dp * (∇(γ) ⋅ nᵥ))dΩfs +
      ∫(dp * γ * (nfs ⋅ nᵥ))dGfs
    )
    jac_t(t, (dtd, dtη, dtu, dtp), (s, γ, v, q)) = (
      ∫(ρf * (v ⋅ dtu))dΩf -
      ∫(αf * dtd * (v ⋅ nᵥ))dΣs -
      ∫(αf * dtη * (v ⋅ nᵥ))dΣfs +
      ∫(ρf * (dtu ⋅ nᵥ) * s)dΩs +
      ∫(ρf * (dtu ⋅ nᵥ) * γ)dΩfs
    )
    jac_tt(t, (dttd, dttη, dttu, dttp), (s, γ, v, q)) = (
      ∫(ρs * hs * (v ⋅ nᵥ) * dttd)dΣs +
      ∫(ρs * hs * s * dttd)dΣs
    )
    l(t, (s, γ, v, q)) = (
      ∫((v ⋅ nᵥ) * -ρf * g + q * 0.0)dΩf +
      ∫(s * -ρf * g)dΩs +
      ∫(γ * -ρf * g)dΩfs -
      ∫(pres_field(t) * v ⋅ (-nᵤ))dΣinlet -
      ∫(pres_field(t) * v ⋅ nᵤ)dΣoutlet
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
    sol_param_fluid_dir = datadir("wsi_2d", "solver_parameters_fluid.xml")
    fluid_block = TrilinosSolve(sol_param_fluid_dir)

    sol_param_solid_dir = datadir("wsi_2d", "solver_parameters_solid.xml")
    solid_block = TrilinosSolve(sol_param_solid_dir)

    fs_block = TrilinosSolve(sol_param_solid_dir)

    coeffs = [
      1.0 1.0 1.0
      0.0 1.0 1.0
      0.0 0.0 1.0
    ]

    bblocks = [
      LinearSystemBlock() LinearSystemBlock() LinearSystemBlock()
      LinearSystemBlock() LinearSystemBlock() LinearSystemBlock()
      LinearSystemBlock() LinearSystemBlock() LinearSystemBlock()
    ]

    prec = BlockTriangularSolver(
      bblocks,
      [solid_block, fs_block, fluid_block],
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
      atol = 1e-6,
      rtol = 1.0e-5,
      verbose = i_am_main(ranks),
    )
    sys_solver = DiscreteDampingSolver(solver, alpha, x_base)

    # ODE Solver
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(sys_solver, dt, ρ∞)

    # Imposing Intitial conditions
    x_initial = interpolate_everywhere(
      [float_solid(t0), free_surface_field(t0), u_field(t0), pres_field(t0)],
      X(t0),
    )
    v_initial = interpolate_everywhere(
      [
        sol_dis_vel_acc(t0),
        fs_ele_vel_acc(t0),
        velocity_vel_acc(t0),
        pressure_vel_acc(t0),
      ],
      X(t0),
    )
    a_initial = interpolate_everywhere(
      [
        sol_dis_vel_acc(t0),
        fs_ele_vel_acc(t0),
        velocity_vel_acc(t0),
        pressure_vel_acc(t0),
      ],
      X(t0),
    )

    # Initializing the solving process
    xₜ = solve(solver_ode, op, t0, tF, (x_initial, v_initial, a_initial))

    @timeit to "Postprocessing Initial Values" begin
      @vtk postprocess_dir = mkpath(datadir("wsi_2d", case, "postprocess"))
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
      @vtk pvd_Σfs = createpvd(
        Σfs,
        ranks,
        "$(postprocess_dir)/freesurface_in_twoD",
      )

      # Storing initial conditions in the pvd file
      @vtk pvd_Ωf[0] = createvtk(
        Ωf,
        "$(vtk_dir)/results_fluids_in_twoD_0",
        cellfields = ["u" => x_initial[3], "pressure" => x_initial[4]],
      )
      @vtk pvd_Σs[0] = createvtk(
        Σs,
        "$(vtk_dir)/results_solids_in_twoD_0",
        cellfields = ["displacement" => x_initial[1]],
      )
      @vtk pvd_Σfs[0] = createvtk(
        Σfs,
        "$(vtk_dir)/results_freesurface_in_twoD_0",
        cellfields = ["FreeSurface" => x_initial[2]],
      )
    end # Ending the postprocessing initial values timer

  end # Ending the model and setup timer

  timestep = 0
  @timeit to "Time Stepping Loop Aggregate" begin
    solve_timer = begin_timed_section!(to, "First solve")
    # Solving and storing results in the pvd file
    with_logger(SimpleLogger(stderr, Logging.Error)) do
      for (tn, (dhn, γhn, uhn, phn)) in xₜ
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
        @vtk pvd_Σfs[tn] = createvtk(
          Σfs,
          "$(vtk_dir)/results_freesurface_in_twoD_$tn",
          cellfields = ["FreeSurface" => γhn],
        )
        end_timed_section!(to, solve_timer)
        i_am_main(ranks) && println("$(timestep + 1) timestep solved")
        timestep += 1
        solve_timer = begin_timed_section!(to, "Subsequent Solves")
      end
    end
    end_timed_section!(to, solve_timer) # Ending the last time step timer
  end # Ending the time stepping loop timer

  @timeit to "Saving pvds" begin
    # Saving the pvd files
    @vtk savepvd(pvd_Ωf)
    @vtk savepvd(pvd_Σs)
    @vtk savepvd(pvd_Σfs)
  end # Ending the saving results timer

  # Displaying the timer summary
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

  GC.gc(true)
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
    freesurface = (
      num_iters = fs_block.log.num_iters,
      residual = fs_block.log.residual,
      solve_time = fs_block.log.solve_time,
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

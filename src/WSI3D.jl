@with_kw struct WSI3D_params
  # Number of MPI processes
  nprocs::Int = 1
  rank::Int = 1

  # Case name
  case::String = "test"

  # Reference dimensional state
  Lref::Float64 = 1.0
  Tref::Float64 = 1.0

  # Physical Parameters
  M::Float64
  τ::Float64

  # Geometrical Parameters (non-dimensionalized)
  H::Float64
  Lm::Float64
  Lf::Float64
  Ly::Float64
  hs::Float64
  meshpath::String = datadir("wsi_3d", "model", "mesh_wsi_3d_1.msh")

  # Damping Parameters (non-dimensionalized)
  Lfd::Float64
  Lfd1::Float64
  Ld::Float64
  Ld1::Float64

  # Temporal parameters (non-dimensionalized)
  ρ∞::Float64
  dt::Float64
  t0::Float64
  tF::Float64

  # Wave Parameters (non-dimensionalized)
  kλ::Float64
  η₀::Float64
  ϕ::Float64
  ω::Float64

  # Robin parameter (non-dimensionalized)
  αf::Float64 = (
    M * (1 - ((2 * ρ∞ - 1) / (1 + ρ∞)))
  ) / (
    dt *
    (1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞)) *
    (1 - (ρ∞ / (1 + ρ∞)))
  ) + 2 * dt / ((1 + ρ∞) * (3 - ρ∞))
  αs::Float64 = 0.0

  # Postprocessing Parameters
  vtkoutput::Bool = false
end

#################################
# Starting FSI Problem Definition
function wsi3d(distribute, parts, params::WSI3D_params)

  # Initializing the timer
  to = TimerOutput()

  @timeit to "Model & Setup" begin

    # ranks
    ranks = distribute(LinearIndices((prod(parts),)))

    # Unpack the parameters
    @unpack_WSI3D_params params

    # Defining the model
    model = UnstructuredDiscreteModel(
      GmshDiscreteModel(ranks, meshpath, renumber = false),
    )

    # Define reference FE
    fe_order = 1
    reffeᵤ = ReferenceFE(raviart_thomas, Float64, fe_order)
    reffeₚ = ReferenceFE(lagrangian, Float64, fe_order)

    # Define triangulations and integration measures
    degree = 2 * fe_order + 2
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf, degree)

    Σs = BoundaryTriangulation(model, tags = ["FloatingSolid"])
    dΣs = Measure(Σs, degree)

    Λ = SkeletonTriangulation(Σs)
    dΛ = Measure(Λ, degree)
    nΛ = get_normal_vector(Λ)

    Σfs = BoundaryTriangulation(model, tags = ["FreeSurface"])
    dΣfs = Measure(Σfs, degree)

    Σinlet = BoundaryTriangulation(model, tags = ["Inlet"])
    dΣinlet = Measure(Σinlet, degree)

    Σoutlet = BoundaryTriangulation(model, tags = ["Outlet"])
    dΣoutlet = Measure(Σoutlet, degree)

    # Define Test FE Functions
    mfs = BlockMultiFieldStyle(3, (1, 1, 3), (1, 2, 3, 4, 5))
    # Free-slip walls: RT/HDiv Dirichlet DOFs constrain only u⋅n.
    # Tangential velocity remains free on both Bed and SideWalls.
    V = TestFESpace(
      model,
      reffeᵤ,
      conformity = :HDiv,
      dirichlet_tags = ["Bed", "SideWalls"],
    )
    Q = TestFESpace(Ωf, reffeₚ, conformity = :L2) # Pressure
    Z = TestFESpace(Σfs, reffeₚ, conformity = :L2) # Free surface elevation
    E = TestFESpace(Σs, reffeₚ, conformity = :H1) # Membrane displacement
    L = TestFESpace(Σs, reffeₚ, conformity = :L2) # Membrane Lagrange multiplier
    Y = MultiFieldFESpace([E, Z, L, V, Q]; style = mfs) # Test Multifield

    # Defining variable fields
    function pres_field(t)
      function inner_function(x)
        if x[1] <= Lfd
          ans = (
            (H - x[3]) -
            ((η₀ * ω * ω) / kλ) *
            (cosh(kλ * x[3]) / sinh(kλ * H)) *
            cos(kλ * x[1] - ω * t + ϕ) *
            ramp(t, Tref)
          )
        else
          ans = (H - x[3])
        end
        return ans
      end
      return inner_function
    end

    function u_field(t)
      ux(t, x) = (
        -ω *
        η₀ *
        (cosh(kλ * x[3]) / sinh(kλ * H)) *
        cos(kλ * x[1] - ω * t + ϕ) *
        ramp(t, Tref)
      )
      uz(t, x) = (
        -ω *
        η₀ *
        (sinh(kλ * x[3]) / sinh(kλ * H)) *
        sin(kλ * x[1] - ω * t + ϕ) *
        ramp(t, Tref)
      )
      function inner_function(x)
        if x[1] <= Lfd
          ans = VectorValue(ux(t, x), 0.0, uz(t, x))
        else
          ans = VectorValue(0.0, 0.0, 0.0)
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
            ((η₀ * ω * ω) / (kλ)) *
            (cosh(kλ * x[3]) / sinh(kλ * H)) *
            cos(kλ * x[1] - ω * t + ϕ) *
            ramp(t, Tref)
          )
        else
          ans = 0.0
        end
        return ans
      end
      return inner_function
    end

    zero_scalar(t) = x -> 0.0
    zero_vector(t) = x -> VectorValue(0.0, 0.0, 0.0)

    # Define Trial FE Functions
    # Zero normal flux on Bed and SideWalls (free slip).
    U = TransientTrialFESpace(V, [zero_vector, zero_vector])
    P = TransientTrialFESpace(Q)
    Zt = TransientTrialFESpace(Z)
    Et = TransientTrialFESpace(E)
    Lt = TransientTrialFESpace(L)
    X = MultiFieldFESpace([Et, Zt, Lt, U, P]; style = mfs)

    # Damping function
    function α(x)
      c1, c2 = 0.1, 20.0
      frac_x = if x[1] >= Ld
        (x[1] - Ld) / (Lf - Ld)
      elseif x[1] < Lfd
        (Lfd - x[1]) / Lfd
      else
        0.0
      end
      return (1.0 - c1 * frac_x^2) *
        (1.0 - expm1(c2 * frac_x^2) / expm1(c2))
    end
    σmax = -log(1e-4) / dt
    σ(x) = σmax * (1.0 - α(x))

    nᵥ = VectorValue(0.0, 0.0, 1.0) # Vertical z direction
    nᵤ = VectorValue(1.0, 0.0, 0.0) # Horizontal x direction

    # Weak form definition
    jac(t, (dη, dζ, dλ, du, dp), (s, γ, μ, v, q)) = (
      # Membrane solve (s): dynamics
      ∫(τ * (∇(dη) ⋅ ∇(s)))dΣs -
      ∫(dλ * s)dΣs +
      # Free surface solve (γ): kinematics
      ∫((-du ⋅ nᵥ) * γ)dΣfs +
      # Fluid solve (μ): Robin interface (Dirichlet)
      ∫(αf * (μ * (du ⋅ nᵥ)))dΣs +
      # Fluid solve (μ): Robin interface (Neumann)
      ∫(τ * (∇(dη) ⋅ ∇(μ)))dΣs -
      ∫(τ * (jump(μ*nΛ)) ⋅ mean(∇(dη)))dΛ -
      ∫(dλ * μ)dΣs +
      # Fluid solve (v):
      ∫(-(∇ ⋅ v) * dp)dΩf +
      ∫(dζ * (v ⋅ nᵥ))dΣfs +
      ∫(dη * (v ⋅ nᵥ))dΣs + ∫(dλ * (v ⋅ nᵥ))dΣs +
      # Fluid solve (q): incompressibility
      ∫((∇ ⋅ du) * q)dΩf +
      # Fluid velocity damping
      ∫(σ * (du ⋅ v))dΩf +
      # Free surface damping
      ∫(σ * (dζ * γ))dΣfs
    )
    jac_t(t, (dtη, dtζ, dtλ, dtu, dtp), (s, γ, μ, v, q)) = (
      # Free surface solve (γ): kinematics
      ∫(dtζ * γ)dΣfs +
      # Fluid solve (μ): Robin interface (Dirichlet)
      ∫(-αf * (dtη * μ))dΣs +
      # Fluid solve (v):
      ∫(v ⋅ dtu)dΩf
    )
    jac_tt(t, (dttη, dttζ, dttλ, dttu, dttp), (s, γ, μ, v, q)) = (
      # Membrane solve (s): inertia
      ∫(M * dttη * s)dΣs +
      # Fluid solve (μ): Robin interface (Neumann)
      ∫(M * dttη * μ)dΣs
    )
    l(t, (s, γ, μ, v, q)) = (
      # Fluid solve (v): gravity and prescribed inlet/outlet pressure
      ∫((v ⋅ nᵥ) * -1.0)dΩf -
      ∫(pres_field(t) * (v ⋅ (-nᵤ)))dΣinlet -
      ∫(pres_field(t) * (v ⋅ nᵤ))dΣoutlet +
      # Fluid velocity damping
      ∫(σ * (u_field(t) ⋅ v))dΩf +
      # Free surface damping
      ∫(σ * (free_surface_field(t) * γ))dΣfs +
      # Zero terms for the other test functions
      ∫(0.0 * s)dΣs + ∫(0.0 * μ)dΣs + ∫(0.0 * q)dΩf + ∫(0.0 * jump(μ))dΛ
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
    sol_param_fluid_dir = datadir("wsi_3d", "solver_parameters_fluid.xml")
    fluid_block = TrilinosSolve(sol_param_fluid_dir)

    sol_param_solid_dir = datadir("wsi_3d", "solver_parameters_solid.xml")
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
      atol = 1e-16,
      rtol = 1.0e-8,
      verbose = i_am_main(ranks),
    )
    sys_solver = SolverWrapper(solver)

    # ODE Solver
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(sys_solver, dt, ρ∞)

    # Initial conditions
    x_initial = interpolate_everywhere(
      [zero_scalar(t0), free_surface_field(t0), zero_scalar(t0), u_field(t0), pres_field(t0)],
      X(t0),
    )
    v_initial = interpolate_everywhere(
      [zero_scalar(t0), zero_scalar(t0), zero_scalar(t0), zero_vector(t0), zero_scalar(t0)],
      X(t0),
    )
    a_initial = interpolate_everywhere(
      [zero_scalar(t0), zero_scalar(t0), zero_scalar(t0), zero_vector(t0), zero_scalar(t0)],
      X(t0),
    )

    # Visualization only: RT velocities are tangentially discontinuous between cells, which
    # VTK shows as speckle at shared vertices. Also output the L2 projection onto continuous P_k.
    Vvis = FESpace(Ωf, ReferenceFE(lagrangian, VectorValue{3,Float64}, fe_order))
    Mvis = assemble_matrix((u, v) -> ∫(u ⋅ v)dΩf, Vvis, Vvis)
    vis_solver = CGSolver(JacobiLinearSolver(); rtol = 1.0e-10, atol = 1.0e-12)
    vis_ns = numerical_setup(symbolic_setup(vis_solver, Mvis), Mvis)
    function project_vis(uh)
      rhs = assemble_vector(v -> ∫(v ⋅ uh)dΩf, Vvis)
      x_vis = allocate_in_domain(Mvis)
      fill!(x_vis, 0)
      solve!(x_vis, vis_ns, rhs)
      return FEFunction(Vvis, x_vis)
    end

    # Initializing the solving process (time steps are solved during iteration).
    xₜ = solve(solver_ode, op, t0, tF, (x_initial, v_initial, a_initial))

    @vtk postprocess_dir = mkpath(datadir("wsi_3d", case, "postprocess"))
    @vtk vtk_dir = mkpath(joinpath(postprocess_dir, "tmp"))

    # Creating pvd files
    @vtk pvd_Ωf = createpvd(Ωf, ranks, joinpath(postprocess_dir, "fluids_in_threeD"))
    @vtk pvd_Σs = createpvd(Σs, ranks, joinpath(postprocess_dir, "solids_in_threeD"))
    @vtk pvd_Σfs = createpvd(Σfs, ranks, joinpath(postprocess_dir, "freesurface_in_threeD"))

    # Storing initial conditions in the pvd file
    @vtk pvd_Ωf[t0] = createvtk(
      Ωf,
      joinpath(vtk_dir, "results_fluids_in_threeD_0"),
      cellfields = [
        "u" => x_initial[4],
        "u_vis" => project_vis(x_initial[4]),
        "pressure" => x_initial[5],
      ],
    )
    @vtk pvd_Σs[t0] = createvtk(
      Σs,
      joinpath(vtk_dir, "results_solids_in_threeD_0"),
      cellfields = ["displacement" => x_initial[1], "membrane_pressure" => x_initial[3]],
    )
    @vtk pvd_Σfs[t0] = createvtk(
      Σfs,
      joinpath(vtk_dir, "results_freesurface_in_threeD_0"),
      cellfields = ["FreeSurface" => x_initial[2]],
    )

  end # Ending the model and setup timer

  timestep = 0
  @timeit to "Time Stepping Loop Aggregate" begin
    solve_timer = begin_timed_section!(to, "First solve")
    with_logger(SimpleLogger(stderr, Logging.Error)) do
      for (tn, state) in xₜ
        ηhn, ζhn, λhn, uhn, phn = state
        @vtk pvd_Ωf[tn] = createvtk(
          Ωf,
          joinpath(vtk_dir, "results_fluids_in_threeD_$tn"),
          cellfields = ["u" => uhn, "u_vis" => project_vis(uhn), "pressure" => phn],
        )
        @vtk pvd_Σs[tn] = createvtk(
          Σs,
          joinpath(vtk_dir, "results_solids_in_threeD_$tn"),
          cellfields = ["displacement" => ηhn, "membrane_pressure" => λhn],
        )
        @vtk pvd_Σfs[tn] = createvtk(
          Σfs,
          joinpath(vtk_dir, "results_freesurface_in_threeD_$tn"),
          cellfields = ["FreeSurface" => ζhn],
        )
        end_timed_section!(to, solve_timer)
        i_am_main(ranks) && println("$(timestep + 1) timestep solved")
        timestep += 1
        solve_timer = begin_timed_section!(to, "Subsequent Solves")
      end
    end
    end_timed_section!(to, solve_timer)
  end

  @timeit to "Saving pvds" begin
    @vtk savepvd(pvd_Ωf)
    @vtk savepvd(pvd_Σs)
    @vtk savepvd(pvd_Σfs)
  end

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
  GC.gc(true)
  return solver_stats
end
# Ending FSI Problem Definition
###############################

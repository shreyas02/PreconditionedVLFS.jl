@with_kw struct Periodic3D_params
    # Number of MPI processes
    nprocs::Int = 1
    rank::Int = 1

    # Case number
    case::String = "test"
    iter::Int = 0

    # Geometric parameters
    H::Float64 = 1.1 # Height of the domain
    Lf::Float64 = 10*π # Length of the entire domain
    Ly::Float64 = 2 # Width of the domain
    hs::Float64 = 0.01 # Floating structure thickness

    # Wave Parameters
    kλ::Float64 = 3.0 # Wave number
    η₀::Float64 = 0.01 # surface elevation

    # Time Numerics
    ρ∞::Float64 = 0.5
    dt::Float64 = 0.1
    t0::Float64 = 0.0
    tF::Float64 = 6.0

    # Physical Parameters
    ρf::Float64 = 1000.0 # Fluid density
    ρs::Float64 = 100.0 # Solid density
    g::Float64 = 9.81 # Acceleration due to gravity

    # Robin parameter
    αf::Float64 =
        (
            ρs*hs*(1-((2 * ρ∞ - 1) / (1 + ρ∞)))
        )/(dt*(1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))*(1-(ρ∞ / (1 + ρ∞))))
    αs::Float64 = 0.0

    # Postprocessing Parameters
    vtkoutput::Bool = false
end

function periodic3D(distribute, parts, params::Periodic3D_params)
    # Initializing the timer
    to = TimerOutput()

    @timeit to "Model & Setup" begin

        # ranks 
        ranks = distribute(LinearIndices((prod(parts),)))

        # Unpack the parameters
        @unpack_Periodic3D_params params

        # Wave parameters
        ω = sqrt(g*kλ*tanh(kλ*H)) # Wave frequency in radians
        ϕ = 0 # wave phase difference

        # Derived parameters
        T = ((ω*ω)/(kλ*kλ))*((ρf)/(kλ*tanh(kλ*H)) + ρs*hs) # initial tension on the floating structure

        # Outer box Coordinates
        domain = (0.0, Lf, 0.0, Ly, 0.0, H)
        nx = trunc(5*Lf);
        ny = trunc(3*Ly);
        nz = 5;
        partition = (nx, ny, nz)
        function f_z(x)
            if x == H
                return H
            end
            i = x / (H/nz)
            return H-H/(1.2^i)
        end
        map_cord(x) = VectorValue(x[1], x[2], f_z(x[3]))
        model = simplexify(
            CartesianDiscreteModel(
                ranks,
                parts,
                domain,
                partition,
                map = map_cord,
                isperiodic = (true, false, false),
            ),
        )
        labels = get_face_labeling(model)
        add_tag_from_tags!(labels, "Bed", [21, 13, 14])
        add_tag_from_tags!(labels, "BedInt", [1, 2, 3, 4, 9, 10])
        add_tag_from_tags!(labels, "Inlet", [25]) # Remove this tag later
        add_tag_from_tags!(labels, "Outlet", [26]) # Remove this tag later
        add_tag_from_tags!(labels, "SideWalls", [11, 12, 23, 24, 17, 18, 19, 20])
        add_tag_from_tags!(labels, "LeftPoint", [15, 5, 7]) # Remove this tag later
        add_tag_from_tags!(labels, "RightPoint", [16, 6, 8]) # Remove this tag later
        add_tag_from_tags!(labels, "FloatingSolid", [22, 11, 12])
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
        FEorder = 2
        reffeᵤ = ReferenceFE(lagrangian, VectorValue{3,Float64}, FEorder)
        reffeₚ = ReferenceFE(lagrangian, Float64, FEorder-1)
        reffeₛ = ReferenceFE(lagrangian, Float64, FEorder)

        # Define triangulation and integration measure
        degree = 2*FEorder + 1

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
            dirichlet_tags = ["Bed", "BedInt", "SideWalls"],
            dirichlet_masks = [
                (false, false, true),
                (false, true, true),
                (false, true, false),
            ],
        )
        Q = TestFESpace(Ωf, reffeₚ, conformity = :H1)
        S = TestFESpace(Ωs, reffeₛ, dirichlet_tags = ["!FloatingSolid"])
        Y = MultiFieldFESpace([S, V, Q]; style = mfs)

        # Define pressure boundary condition functions 
        ux(t, x) = ω*η₀*(cosh(kλ*(x[3]))/(sinh(kλ*H)))*cos(kλ*x[1] - ω*t + ϕ)
        uz(t, x) = ω*η₀*(sinh(kλ*(x[3]))/(sinh(kλ*H)))*sin(kλ*x[1] - ω*t + ϕ)
        u_field(t) = x -> VectorValue(ux(t, x), 0.0, uz(t, x))
        pres_field(t) =
            x->ρf*g*(H - x[3]) +
               ρf*((η₀*ω*ω)/(kλ))*(cosh(kλ*(x[3]))/(sinh(kλ*H)))*cos(kλ*x[1] - ω*t + ϕ)
        float_solid(t) = x->η₀*cos(kλ*x[1] - ω*t + ϕ)
        zero_scalar(t) = x -> 0.0

        # Dual use (vel/acc) terms
        pressure_VelAcc(t) = x -> 0.0 # Initial vel/acc terms for pressure
        velocity_VelAcc(t) = x -> VectorValue(0.0, 0.0, 0.0) # Initial vel/acc terms for velocity
        sol_dis_VelAcc(t) = x -> 0.0 # Initial vel/acc terms for solid displacement

        # Define Trial FE Functions
        U = TransientTrialFESpace(V, [u_field, u_field, u_field]) # Trial function for fluid velocity
        P = TransientTrialFESpace(Q) # Trial function for Pressure
        D = TransientTrialFESpace(S, [zero_scalar]) # Trial function for membrane Displacement
        X = MultiFieldFESpace([D, U, P]; style = mfs) # Trial Multifield

        # Damping function
        α(x) = 1.0
        αvec(x) = VectorValue(α(x), α(x), α(x))
        alpha = interpolate_everywhere([α, αvec, α], X(0.0))
        x_base(t) =
            interpolate_everywhere([float_solid(t), u_field(t), pres_field(t)], X(t))

        nᵥ=VectorValue(0.0, 0.0, 1.0) # Along the Z direction
        nᵤ=VectorValue(1.0, 0.0, 0.0) # ALong the X direction

        # Weak form definition
        jac(t, (dd, du, dp), (s, v, q)) =
            ∫(-(∇⋅v)*dp)dΩf +
            ∫((∇⋅du)*q)dΩf + # Fluid problem (Interior)
            ∫(T*(((∇(dd)⋅nᵤ)*(∇(v⋅nᵥ)⋅nᵤ))))dΣs + # Fluid problem (Neumann) 
            ∫(αf*(du⋅nᵥ)*(v⋅nᵥ))dΣs + # Fluid problem (Dirichlet)
            ∫(T*(∇(dd)⋅nᵤ)*(∇(s)⋅nᵤ))dΣs - ∫(dp*(∇(s)⋅nᵥ))dΩs  # Solid problem 
        jac_t(t, (dtd, dtu, dtp), (s, v, q)) =
            ∫(ρf*(v⋅dtu))dΩf - # Fluid problem (Interior)
            ∫(αf*(dtd)*(v⋅nᵥ))dΣs + # Fluid problem (Dirichlet)
            ∫(ρf*(dtu⋅nᵥ)*s)dΩs  # Solid problem
        jac_tt(t, (dttd, dttu, dttp), (s, v, q)) =
            ∫(ρs*hs*(v⋅nᵥ)*dttd)dΣs + # Fluid problem (Neumann)
            ∫(ρs*hs*s*dttd)dΣs # Solid problem
        l(t, (s, v, q)) =
            ∫((v⋅nᵥ)*-ρf*g + (q*0.0))dΩf + # Fluid problem
            ∫(s*-ρf*g)dΩs # Solid problem

        # Build affine FE operator
        op = TransientLinearFEOperator((jac, jac_t, jac_tt), l, X, Y)

        # System Solver Definition
        sol_param_fluid_dir = datadir("periodic3d", "solver_parameters_fluid.xml")
        fluidBlock = TrilinosSolve(sol_param_fluid_dir, 1000)

        sol_param_solid_dir = datadir("periodic3d", "solver_parameters_solid.xml")
        solidBlock = TrilinosSolve(sol_param_solid_dir, 1000)

        coeffs = [
            1.0 1.0;
            0.0 1.0;
        ]

        bblocks = [
            LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock();
        ]

        P = BlockTriangularSolver(
            bblocks,
            [solidBlock, fluidBlock],
            coeffs,
            :upper,
        )
        solver = FGMRESSolver(
            20,
            P;
            Pl = nothing,
            restart = false,
            m_add = 1,
            maxiter = 1000,
            atol = 1e-6,
            rtol = 1.e-5,
            verbose = i_am_main(ranks),
        )
        SysSolver = DiscreteDampingSolver(solver, alpha, x_base)

        # ODE Solver
        solver_ode = Gridap.ODEs.GeneralizedAlpha2(SysSolver, dt, ρ∞);

        # Imposing Intitial conditions
        x_initial =
            interpolate_everywhere([float_solid(t0), u_field(t0), pres_field(t0)], X(t0))
        v_initial = interpolate_everywhere(
            [sol_dis_VelAcc(t0), velocity_VelAcc(t0), pressure_VelAcc(t0)],
            X(t0),
        )
        a_initial = interpolate_everywhere(
            [sol_dis_VelAcc(t0), velocity_VelAcc(t0), pressure_VelAcc(t0)],
            X(t0),
        )

        # Initializing the solving process
        xₜ = solve(solver_ode, op, t0, tF, (x_initial, v_initial, a_initial))

        @timeit to "Postprocessing Initial Values" begin
            @vtk postprocess_dir = mkpath("$(datadir("periodic3d", string(case), "postprocess"))")
            @vtk vtk_dir = mkpath("$(postprocess_dir)/tmp")

            # Creating pvd files
            @vtk pvd_Ωf = createpvd(Ωf, ranks, "$(postprocess_dir)/fluids_in_threeD")
            @vtk pvd_Σs = createpvd(Σs, ranks, "$(postprocess_dir)/solids_in_threeD")

            # Storing initial conditions in the pvd file 
            @vtk pvd_Ωf[0] = createvtk(
                Ωf,
                "$(vtk_dir)/results_fluids_in_threeD_0",
                cellfields = ["u" => x_initial[2], "pressure" => x_initial[3]],
            )
            @vtk pvd_Σs[0] = createvtk(
                Σs,
                "$(vtk_dir)/results_solids_in_threeD_0",
                cellfields = ["displacement" => x_initial[1]],
            )
        end # Ending the postprocessing initial values timer

    end # Ending the model and setup timer

    timestep = 0;
    @timeit to "Time Stepping Loop Aggregate" begin
        solve_timer = begin_timed_section!(to, "First solve")
        # Solving and storing results in the pvd file 
        with_logger(SimpleLogger(stderr, Logging.Error)) do  # Warnings in the vtu files. Need to fix later 
            for (tn, (dhn, uhn, phn)) in xₜ
                @vtk pvd_Ωf[tn] = createvtk(
                    Ωf,
                    "$(vtk_dir)/results_fluids_in_threeD_$tn",
                    cellfields = ["u" => uhn, "pressure" => phn],
                )
                @vtk pvd_Σs[tn] = createvtk(
                    Σs,
                    "$(vtk_dir)/results_solids_in_threeD_$tn",
                    cellfields = ["displacement" => dhn],
                )
                end_timed_section!(to, solve_timer)
                i_am_main(ranks) && println("$(timestep+1) timestep solved")
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
    outer_residuals = outer_residuals[outer_iter_array.+1] # Adjusting for 1-based indexing

    return (
        fluid = (num_iters = fluidBlock.log.num_iters(),
                 residual = fluidBlock.log.residual(),
                 solve_time = fluidBlock.log.solve_time()),
        solid = (num_iters = solidBlock.log.num_iters(),
                 residual = solidBlock.log.residual(),
                 solve_time = solidBlock.log.solve_time()),
        outer = (iter_array = outer_iter_array,
                residuals = outer_residuals,
                timer = to))
end
# Ending FSI Problem Definition 
###############################
@with_kw struct WSI3D_params
    # Number of MPI processes
    nprocs::Int = 1
    rank::Int = 1

    # Case name
    case::String = "test"

    # Geometrical Parameters
    H::Float64
    Lm::Float64
    Lf::Float64
    Ly::Float64
    hs::Float64
    meshpath::String = datadir("wsi_3d", "model", "mesh_wsi_3d_1.msh")
    
    # Damping Parameters
    Lfd::Float64
    Lfd1::Float64
    Ld::Float64
    Ld1::Float64

    # Temporal parameters
    ρ∞::Float64
    dt::Float64
    t0::Float64
    tF::Float64

    # Physical Parameters
    ρf::Float64
    ρs::Float64
    g::Float64
    T::Float64

    # Wave Parameters
    kλ::Float64
    η₀::Float64
    ϕ::Float64
    ω::Float64

    # Robin parameter
    αf::Float64 =
        (
            ρs*hs*(1-((2 * ρ∞ - 1) / (1 + ρ∞)))
        )/(dt*(1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))*(1-(ρ∞ / (1 + ρ∞))))
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
    model = UnstructuredDiscreteModel(GmshDiscreteModel(ranks, meshpath, renumber=false))
    labels = get_face_labeling(model)
    Geometry.add_tag_from_tags_complementary!(labels,"!!FreeSurface",["FreeSurface"])
    Geometry.add_tag_from_tags_setdiff!(labels,"!FreeSurface",["!!FreeSurface"],["LeftPoint","RightPoint"])
    Geometry.add_tag_from_tags_complementary!(labels,"!FloatingSolid",["FloatingSolid"])


    # Define reference FE (Q2/P1(disc) pair)
    FEorder = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{3,Float64}, FEorder)
    reffeₚ = ReferenceFE(lagrangian, Float64, FEorder-1)
    reffeₛ = ReferenceFE(lagrangian, Float64, FEorder)
    
    # Define triangulation and integration measure
    degree = 2*FEorder + 1
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf, degree)

    Ωfs, Ωs = get_triangulations(model,"FreeSurface","FloatingSolid") 
    
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
    mfs = BlockMultiFieldStyle(3,(1,1,2),(1,2,3,4))
    V = TestFESpace(Ωf, reffeᵤ, dirichlet_tags = ["Bed","BedInt","SideWalls","LeftPoint", "RightPoint"],
                    dirichlet_masks = [(false, false, true),(false, true, true),(false, true, false),(true, true, true),(true, true, true)]) # Test function for fluid velocity
    Q = TestFESpace(Ωf, reffeₚ, conformity=:H1) # Test function for Pressure 
    S = TestFESpace(Ωs,reffeₛ, dirichlet_tags = ["!FloatingSolid"]) # Test function for Solid Displacement
    Sfs = TestFESpace(Ωfs, reffeₛ, dirichlet_tags = ["LeftPoint", "RightPoint", "!FreeSurface"]) # Test function for free surface elevation
    Y = MultiFieldFESpace([S, Sfs, V, Q]; style = mfs) # Test Multifield

    Gs = select_triangulation(Ωs, Ωfs)
    dGs = Measure(Gs, degree)
    Gfs = select_triangulation(Ωfs, Ωs)
    dGfs = Measure(Gfs, degree)
    ns = get_normal_vector(Gs)
    nfs = get_normal_vector(Gfs)

    ## Defining variable fields -
    function ramp_t(t)
        if t <= 2.0 
            frac = (t - 0.0)/2.0
            ans = 3*frac^2 - 2*frac^3
            return ans
        else
            return 1.0
        end
    end
    function u_field(t)
        ux(t,x) = -ω*η₀*(cosh(kλ*(x[3]))/(sinh(kλ*H)))*cos(kλ*x[1] - ω*t + ϕ) * ramp_t(t)
        uy(t,x) = 0.0
        uz(t,x) = -ω*η₀*(sinh(kλ*(x[3]))/(sinh(kλ*H)))*sin(kλ*x[1] - ω*t + ϕ) * ramp_t(t)
        function inner_function(x)
            if(x[1] <= Lfd)
                ans = VectorValue(ux(t,x), uy(t,x), uz(t,x))
            else
                ans = VectorValue(0.0, 0.0, 0.0)
            end
            return ans
        end
        return inner_function
    end
    function pres_field(t)
        function inner_function(x)
            if(x[1] <= Lfd)
                ans = ρf*g*(H - x[3]) - ((ρf*η₀*ω*ω)/(kλ))*(cosh(kλ*(x[3]))/(sinh(kλ*H)))*cos(kλ*x[1] - ω*t + ϕ) * ramp_t(t)
            else
                ans = ρf*g*(H - x[3])
            end
            return ans
        end
        return inner_function
    end
    function freeSurface_field(t)
        function inner_function(x)
            if(x[1] <= Lfd)
                ans = -1.0*((η₀*ω*ω)/(kλ*g))*(cosh(kλ*(x[3]))/(sinh(kλ*H)))*cos(kλ*x[1] - ω*t + ϕ) * ramp_t(t)
            else
                ans = 0.0
            end
            return ans
        end
        return inner_function
    end
    float_solid(t) = x -> 0.0 # Explicit boundary/initial condition for solid displacement
  
    # Dual use (vel/acc) terms
    pressure_VelAcc(t) = x -> 0.0 # Initial vel/acc terms for pressure
    velocity_VelAcc(t) = x -> VectorValue(0.0, 0.0, 0.0) # Initial vel/acc terms for velocity
    sol_dis_VelAcc(t) = x -> 0.0 # Initial vel/acc terms for solid displacement 
    fs_ele_VelAcc(t) = x -> 0.0 # Initial vel/acc terms for free surface elevation

    # Define Trial FE Functions
    U = TransientTrialFESpace(V, [velocity_VelAcc,velocity_VelAcc,velocity_VelAcc,u_field,u_field]) # Trial function for fluid velocity
    P = TransientTrialFESpace(Q) # Trial function for Pressure
    D = TransientTrialFESpace(S, [sol_dis_VelAcc]) # Trial function for Solid Displacement
    Dfs = TransientTrialFESpace(Sfs, [freeSurface_field, freeSurface_field, fs_ele_VelAcc]) # Trial function for free surface elevation
    X = MultiFieldFESpace([D, Dfs, U, P]; style = mfs) # Trial Multifield

    # Time stepping parameters for generalized alpha method

    # Damping function
    function α(x)
        c1 = 0.1
        c2 = 20.0
        function frac(x,Ld,Ld1, Lfd, Lfd1)
            if(x[1] ≈ Ld1 || x[1] > Ld1)
                return  1.0 # 100% outlet damping after Ld1
            elseif(x[1] >= Ld)
                return (x[1] - Ld)/(Ld1 - Ld) # Outlet damping between Ld and Ld1
            elseif(x[1] ≈ Lfd1 || x[1] < Lfd1)
                return  1.0 # 100% inlet damping before Lfd1
            elseif(x[1] <= Lfd)
                return (Lfd - x[1])/(Lfd - Lfd1) # Inlet damping between Lfd1 and Lfd
            else
                return 0.0 # No damping in the middle region
            end
        end
        return (1.0 - c1*frac(x,Ld,Ld1,Lfd,Lfd1)*frac(x,Ld,Ld1,Lfd,Lfd1))*(1.0 - (1.0 - exp(c2*frac(x,Ld,Ld1,Lfd,Lfd1)*frac(x,Ld,Ld1,Lfd,Lfd1)))/(1.0 - exp(c2)))
    end
    αvec(x) = VectorValue(α(x),α(x),α(x))
    alpha = interpolate_everywhere([α,α,αvec,α], X(0.0))
    x_base(t) = interpolate_everywhere([float_solid(t),freeSurface_field(t),u_field(t),pres_field(t)],X(t))

    nᵥ=VectorValue(0.0,0.0,1.0) # Along the Z direction
    nᵤ=VectorValue(1.0,0.0,0.0) # ALong the X direction

    # Weak form definition
    jac(t,(dd,dη,du,dp),(s,γ,v,q)) = ∫(-(∇⋅v)*dp)dΩf + ∫((∇⋅du)*q)dΩf + # Fluid problem (Interior)
            ∫(T*(((∇(dd)⋅nᵤ)*(∇(v⋅nᵥ)⋅nᵤ))))dΣs + ∫(ρf*g*dη*(v⋅nᵥ))dΣfs + # Fluid problem (Neumann) 
            ∫(αf*(du⋅nᵥ)*(v⋅nᵥ))dΣs + ∫(αf*(du⋅nᵥ)*(v⋅nᵥ))dΣfs + # Fluid problem (Dirichlet)
            ∫(T*(∇(dd)⋅nᵤ)*(∇(s)⋅nᵤ))dΣs - ∫(dp*(∇(s)⋅nᵥ))dΩs + ∫(dp*s*(ns⋅nᵥ))dGs +# Solid problem 
            ∫(ρf*g*dη*γ)dΣfs - ∫(dp*(∇(γ)⋅nᵥ))dΩfs + ∫(dp*γ*(nfs⋅nᵥ))dGfs  # Free surface problem
    jac_t(t,(dtd,dtη,dtu,dtp),(s,γ,v,q)) = ∫(ρf*(v⋅dtu))dΩf - # Fluid problem (Interior)
            ∫(αf*(dtd)*(v⋅nᵥ))dΣs - ∫(αf*dtη*(v⋅nᵥ))dΣfs + # Fluid problem (Dirichlet)
            ∫(ρf*(dtu⋅nᵥ)*s)dΩs + # Solid problem
            ∫(ρf*(dtu⋅nᵥ)*γ)dΩfs # Free surface problem
    jac_tt(t,(dttd,dttη,dttu,dttp),(s,γ,v,q)) = ∫(ρs*hs*(v⋅nᵥ)*dttd)dΣs + # Fluid problem (Neumann)
            ∫(ρs*hs*s*dttd)dΣs # Solid problem
    l(t,(s,γ,v,q)) = ∫((v⋅nᵥ)*-ρf*g + (q*0.0))dΩf + # Fluid problem
            ∫(s*-ρf*g)dΩs + # Solid problem
            ∫(γ*-ρf*g)dΩfs - # Free surface problem
            ∫(pres_field(t)*v⋅(-nᵤ))dΣinlet - ∫(pres_field(t)*v⋅(nᵤ))dΣoutlet # Pressure BCs

    # Build affine FE operator
    op = TransientLinearFEOperator((jac, jac_t, jac_tt), l, X, Y)   

    # System Solver Definition
    sol_param_fluid_dir = datadir("wsi_3d", "solver_parameters_fluid.xml")
    fluidBlock = TrilinosSolve(sol_param_fluid_dir, 1000)

    sol_param_solid_dir = datadir("wsi_3d", "solver_parameters_solid.xml")
    solidBlock = TrilinosSolve(sol_param_solid_dir, 1000)

    fsBlock = TrilinosSolve(sol_param_solid_dir, 1000)

    coeffs = [1.0 1.0 1.0;
              0.0 1.0 1.0;
              0.0 0.0 1.0]

    bblocks = [LinearSystemBlock() LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock() LinearSystemBlock();]

    P = BlockTriangularSolver(bblocks,[solidBlock,fsBlock,fluidBlock],coeffs,:upper)
    solver = FGMRESSolver(20,P;Pl=nothing,restart=false,m_add=1,maxiter=1000,atol=1e-6,rtol=1.e-5,verbose=i_am_main(ranks))
    SysSolver = DiscreteDampingSolver(solver, alpha, x_base)

    # ODE Solver
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(SysSolver, dt, ρ∞)

    # Imposing Intitial conditions
    x_initial = interpolate_everywhere([float_solid(t0), freeSurface_field(t0), u_field(t0), pres_field(t0)], X(t0))
    v_initial = interpolate_everywhere([sol_dis_VelAcc(t0),fs_ele_VelAcc(t0), velocity_VelAcc(t0), pressure_VelAcc(t0)], X(t0))
    a_initial = interpolate_everywhere([sol_dis_VelAcc(t0), fs_ele_VelAcc(t0), velocity_VelAcc(t0), pressure_VelAcc(t0)], X(t0))

    # Initializing the solving process
    xₜ = solve(solver_ode, op, t0, tF, (x_initial,v_initial,a_initial))
    
    @timeit to "Postprocessing Initial Values" begin
    @vtk postprocess_dir = mkpath("$(datadir("wsi_3d",case,"postprocess"))")
    @vtk vtk_dir = mkpath("$(postprocess_dir)/tmp")

    # Creating pvd files
    @vtk pvd_Ωf = createpvd(Ωf,ranks,"$(postprocess_dir)/fluids_in_threeD")
    @vtk pvd_Σs = createpvd(Σs,ranks,"$(postprocess_dir)/solids_in_threeD")
    @vtk pvd_Σfs = createpvd(Σfs,ranks,"$(postprocess_dir)/freesurface_in_threeD")

    # Storing initial conditions in the pvd file 
    @vtk pvd_Ωf[0] = createvtk(Ωf,"$(vtk_dir)/results_fluids_in_threeD_0",cellfields = ["u" => x_initial[3], "pressure" => x_initial[4]],)
    @vtk pvd_Σs[0] = createvtk(Σs,"$(vtk_dir)/results_solids_in_threeD_0",cellfields = ["displacement" => x_initial[1]],)
    @vtk pvd_Σfs[0] = createvtk(Σfs,"$(vtk_dir)/results_freesurface_in_threeD_0" * ".vtu",cellfields = ["FreeSurface" => x_initial[2]],)
    end # Ending the postprocessing initial values timer

    end # Ending the model and setup timer

    timestep = 0;
    @timeit to "Time Stepping Loop Aggregate" begin
    solve_timer = begin_timed_section!(to, "First solve")
    # Solving and storing results in the pvd file 
    with_logger(SimpleLogger(stderr, Logging.Error)) do # Warnings in the vtu files. Need to fix later 
        for (tn, (dhn,γhn,uhn,phn)) in xₜ
            @vtk pvd_Ωf[tn] = createvtk(Ωf, "$(vtk_dir)/results_fluids_in_threeD_$tn", cellfields=["u" => uhn, "pressure" => phn])
            @vtk pvd_Σs[tn] = createvtk(Σs, "$(vtk_dir)/results_solids_in_threeD_$tn", cellfields=["displacement" => dhn])
            @vtk pvd_Σfs[tn] = createvtk(Σfs, "$(vtk_dir)/results_freesurface_in_threeD_$tn", cellfields=["FreeSurface" => γhn])
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
    outer_residuals = outer_residuals[outer_iter_array.+1] # Adjusting for 1-based indexing

    return (
        fluid = (num_iters = fluidBlock.log.num_iters(),
                 residual = fluidBlock.log.residual(),
                 solve_time = fluidBlock.log.solve_time()),
        solid = (num_iters = solidBlock.log.num_iters(),
                 residual = solidBlock.log.residual(),
                 solve_time = solidBlock.log.solve_time()),
        freesurface = (num_iters = fsBlock.log.num_iters(),
                       residual = fsBlock.log.residual(),
                       solve_time = fsBlock.log.solve_time()),
        outer = (iter_array = outer_iter_array,
                residuals = outer_residuals,
                timer = to))
end
# Ending FSI Problem Definition 
###############################

@with_kw struct ToyRichardson_params
    # Geometrical Parameters
    length::Float64 = 6.0 # Length of the domain
    height::Float64 = 1.0 # Height of the domain
    hs::Float64 = 0.1 # Thickness of the solid

    # Physical Parameters
    ρf::Float64 = 1.0 # Fluid Density
    ρs::Float64 = 45.0 # Solid Density
    L::Float64 = 1e5 # L = hs*E/(R²(1-ν²))

    # Time Numerics
    ρ∞::Float64 = 0.5
    dt::Float64 = 0.001
    t0::Float64 = 0.0
    tF::Float64 = 10.0*dt

    # Robin Parameters
    αf::Float64 = (ρs*hs*(1-((2 * ρ∞ - 1) / (1 + ρ∞))))/(dt*(1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))*(1-(ρ∞ / (1 + ρ∞)))) + (L*dt*((1 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))^2 / 4))/(1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))
    αs::Float64 = 0.0

    # Postprocessing Parameters
    vtkoutput::Bool = false
end

# Functions to run the code
function toyrichardson(ω::Float64, params::ToyRichardson_params)

    # Unpack the parameters
    @unpack_ToyRichardson_params params

    # Defining the Mesh
    domain = (0,length,0,height)
    partition = (length*7,height*7)
    model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain, partition))

    # Define tags in the model
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"Membrane",[6,])
    add_tag_from_tags!(labels,"Left",[1,7,3])
    add_tag_from_tags!(labels,"Right",[2,4,8])
    add_tag_from_tags!(labels,"Bottom",[1,5,2])
    add_tag_from_tags!(labels,"BCMembrane",[3,4])
    Geometry.add_tag_from_tags_complementary!(labels,"!!Membrane",["Membrane"])
    Geometry.add_tag_from_tags_setdiff!(labels,"!Membrane",["!!Membrane"],["BCMembrane"])

    # Define reference FE (P2/P1 taylor hood elements for FEorder = 2)
    FEorder = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},FEorder)
    reffeₚ = ReferenceFE(lagrangian,Float64,FEorder-1)
    reffeₛ = ReferenceFE(lagrangian, Float64, FEorder)

    # Define triangulation and integration measure
    degree = 2*FEorder + 1

    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,degree)

    Ωs = Triangulation(model)
    dΩs = Measure(Ωs, degree)
    
    Σ = BoundaryTriangulation(model,tags="Membrane")
    dΣ = Measure(Σ,degree)
    
    Σleft = BoundaryTriangulation(model,tags="Left")
    dΣleft = Measure(Σleft,degree)
    
    Σright = BoundaryTriangulation(model,tags="Right")
    dΣright = Measure(Σright,degree)

    # Define Test FE Functions
    mfs = BlockMultiFieldStyle(2, (2, 1), (1, 2, 3))
    V = TestFESpace(Ωf,reffeᵤ,dirichlet_tags = ["Bottom","BCMembrane"],
                        dirichlet_masks =[(false,true),(false,true)] ,conformity=:H1)
    Q = TestFESpace(Ωf,reffeₚ, conformity=:H1)
    S = TestFESpace(Ωs,reffeₛ,dirichlet_tags = ["!Membrane","BCMembrane"], conformity=:H1)
    Y = MultiFieldFESpace([V,Q,S]; style = mfs)

    # Define trial FESpaces from Dirichlet functions
    gu(t) = x -> VectorValue(0.0,0.0) # Velocity BCs
    gd(t) = x -> 0.0 # Displacement BCs
    gpr(t) = x -> 0.0 # Right Pressure BCs
    function gpl(t) # Left Pressure BCs
        function inner_function(x)
            if t > 0.005
                return 0.0
            else
                return 20000
            end
        end
        return inner_function
    end
    U = TransientTrialFESpace(V,[gu,gu])
    P = TransientTrialFESpace(Q)
    D = TransientTrialFESpace(S, [gd,gd])
    X = TransientMultiFieldFESpace([U,P,D], style = mfs)

    # Initial condition functions
    h_u(t) = x -> VectorValue(0.0,0.0) # Velocity
    h_d(t) = x -> 0.0 # Displacement
    function h_p(t) # Left Pressure BCs
        function inner_function(x)
            if t > 0.005
                return 0.0
            else
                if x[1] == 0
                    return 20000
                else
                    return 0
                end
            end
        end
        return inner_function
    end

    # Initial Velocity and Acceleration terms 
    v_p(t) = x -> 0.0 # Pressure
    v_u(t) = x -> VectorValue(0.0,0.0) # Velocity
    v_d(t) = x -> 0.0 # Displacement

    nᵤ=VectorValue(1.0,0.0)
    nᵥ=VectorValue(0.0,1.0)

    # Weak form definition

    jac(t,(du,dp,dd),(v,q,s)) = ∫(-(∇⋅v)*dp + (∇⋅du)*q)dΩf + # Fluid Problem (Interior)
                    ∫(L*dd*(v⋅nᵥ))dΣ + ∫(αf*((du⋅nᵥ)*(v⋅nᵥ)))dΣ + # Fluid Coupling (Robin)
                    ∫(L*dd*s)dΣ + ∫(-αs*((du⋅nᵥ)*s))dΣ - ∫(dp*(∇(s)⋅nᵥ))dΩs # Solid Coupling (Robin)
    jac_t(t,(dtu,dtp,dtd),(v,q,s)) = ∫(ρf*v⋅dtu)dΩf + # Fluid Problem (Interior)
                    ∫(ρf*s⋅(dtu⋅nᵥ))dΩs -
                    ∫(αf*(dtd*(v⋅nᵥ)))dΣ - # Fluid Coupling (Robin)
                    ∫(-αs*(dtd*s))dΣ # Solid Coupling (Robin)
    jac_tt(t,(dttu,dttp,dttd),(v,q,s)) = ∫(ρs*hs*(v⋅nᵥ)*dttd)dΣ + # Fluid Coupling (Robin)
                    ∫(ρs*hs*s*dttd)dΣ # Solid Coupling (Robin)
    l(t,(v,q,s)) = ∫(v⋅VectorValue(0.0, 0.0))dΩf + ∫(q*0.0)dΩf + ∫(s*0.0)dΣ - # Fluid Problem (Interior)
                    ∫(gpl(t)*v⋅-nᵤ)dΣleft - ∫(gpr(t)*v⋅nᵤ)dΣright # Fluid BCs (Neumann)

    # Build affine FE operator
    op = TransientLinearFEOperator((jac, jac_t, jac_tt), l, X, Y)

    # Creation of the relaxation parameter
    _get_ones_scalar(t) = x -> 1.0
    _get_ones_vector(t) = x -> VectorValue(1.0,1.0)
    _get_relax_scalars(t) = x-> ω
    ω_vec = get_free_dof_values(interpolate_everywhere([_get_ones_vector(t0) , _get_ones_scalar(t0), _get_relax_scalars(t0)], X(t0)))

    # System Solver Definition 
    solver_blocks = LUSolver()
    bblocks = [LinearSystemBlock() LinearSystemBlock();
                LinearSystemBlock() LinearSystemBlock()]
    coeffs = [1.0 0.0;
              1.0 1.0]
    P = BlockTriangularSolver(bblocks,[solver_blocks,solver_blocks],coeffs,:lower)
    solver = RichardsonLinearSolver(ω_vec,51;Pl=P,rtol=1e-6,atol=1e-5,verbose=true)

    # ODE Solver
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(solver,dt,ρ∞);

    # Imposing Intitial conditions
    x_initial = interpolate_everywhere([h_u(t0) , h_p(t0), h_d(t0)], X(t0))
    v_initial = interpolate_everywhere([v_u(t0) , v_p(t0), v_d(t0)], X(t0))
    a_initial = interpolate_everywhere([v_u(t0) , v_p(t0), v_d(t0)], X(t0))

    # Initializing the solving process
    xₜ = solve(solver_ode, op, t0, tF, (x_initial,v_initial,a_initial))

    # Post processing begin
    @vtk postprocess_dir = datadir("toy_richardson", "postprocess")
    @vtk vtk_dir = mkpath("$(postprocess_dir)/tmp")

    # Creating pvd files
    @vtk pvd_Ωf = paraview_collection("$(postprocess_dir)/Fluids_in_2D")
    @vtk pvd_Σ = paraview_collection("$(postprocess_dir)/Solids_in_2D")

    # Storing initial conditions in the pvd file 
    @vtk pvd_Ωf[0] = createvtk(Ωf, "$(vtk_dir)/results_Fluids_in_2D_0", cellfields=["u" => x_initial[1] , "pressure" => x_initial[2]])
    @vtk pvd_Σ[0] = createvtk(Σ, "$(vtk_dir)/results_Solids_in_2D_0", cellfields=["displacement" => x_initial[3]])

    timestep = 0
    for (tn, (uhn,phn,dhn)) in xₜ
        @vtk pvd_Ωf[tn] = createvtk(Ωf, "$(vtk_dir)/results_Fluids_in_2D_$tn",cellfields=["u" => uhn, "pressure" => phn])
        @vtk pvd_Σ[tn] = createvtk(Σ, "$(vtk_dir)/results_Solids_in_2D_$tn", cellfields=["displacement" => dhn])
        println("Solving for time step $timestep")
        timestep += 1
    end
    @vtk vtk_save(pvd_Ωf)
    @vtk vtk_save(pvd_Σ)
   
    # Return convergence and timer data for the solver at time tF
    outer_num_iter = solver.log.num_iters
    outer_iter_array = collect(0:outer_num_iter)
    outer_residuals = solver.log.residuals
    outer_residuals = outer_residuals[outer_iter_array.+1] # Adjusting for 1-based indexing

    return outer_num_iter, outer_iter_array, outer_residuals
end
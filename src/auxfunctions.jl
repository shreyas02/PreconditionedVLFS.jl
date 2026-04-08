###################
##### Solvers #####
###################

#######################
## Begin Trilinos Solve 
struct TrilinosSolve <: Gridap.Algebra.LinearSolver
    parameter_file::String
    max_entries_per_row::Int
    log::Base.RefValue{Union{Nothing,GridapTrilinos.SolverResultAllocated }}
end

function TrilinosSolve(parameter_file::AbstractString, max_entries_per_row::Integer=100)
    return TrilinosSolve(
        String(parameter_file),
        Int(max_entries_per_row),
        Ref{Union{Nothing,GridapTrilinos.SolverResultAllocated}}(nothing),
    )
end

log(solver::TrilinosSolve) = solver.log[]

function Base.getproperty(solver::TrilinosSolve, name::Symbol)
    if name === :log
        value = getfield(solver, :log)[]
        value === nothing && error("Solver log is not available yet. Run the solve first.")
        return value
    end
    return getfield(solver, name)
end

struct TrilinosSolveSymbolicSetup <: Gridap.Algebra.SymbolicSetup
    solver::TrilinosSolve
end

mutable struct TrilinosSolveNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
    solver::TrilinosSolve
    A::T
end

function Gridap.Algebra.symbolic_setup(solver::TrilinosSolve, mat::AbstractMatrix)
    return TrilinosSolveSymbolicSetup(solver)
end

function Gridap.Algebra.numerical_setup(ss::TrilinosSolveSymbolicSetup, A::AbstractMatrix)
    return TrilinosSolveNumericalSetup(ss.solver, A)
end

function Gridap.Algebra.numerical_setup(
    ss::TrilinosSolveSymbolicSetup,
    A::AbstractMatrix,
    x::AbstractVector,
)
    return TrilinosSolveNumericalSetup(ss.solver, A)
end

function Gridap.Algebra.numerical_setup!(ns::TrilinosSolveNumericalSetup, A::AbstractMatrix)
    return ns
end

function Gridap.Algebra.solve!(
    x::AbstractVector,
    ns::TrilinosSolveNumericalSetup,
    b::AbstractVector,
)
    map(
        ns.A.matrix_partition,
        ns.A.row_partition,
        ns.A.col_partition,
        b.vector_partition,
        x.vector_partition,
        x.index_partition,
    ) do LocMat, RowMap, ColMap, LocRhs, LocSoln, xmap
        OwnRowMap = RowMap[own_to_local(RowMap)];
        getfield(ns.solver, :log)[] = GridapTrilinos.TrilinosParallel(
            LocMat.nzval,
            LocMat.rowval .- 1,
            LocMat.colptr .- 1,
            LocRhs,
            LocSoln,
            OwnRowMap .- 1,
            ColMap .- 1,
            size(ns.A)[1],
            size(OwnRowMap)[1],
            size(LocMat)[2],
            own_to_local(xmap) .- 1,
            own_to_local(RowMap) .- 1,
            ns.solver.max_entries_per_row,
            ns.solver.parameter_file,
        )
    end
    return ns
end
## End Trilinos Solve
#####################

################################
# Begin DiscreteDampingSolver setup 
struct DiscreteDampingSolver <: LinearSolver
    ls::LinearSolver
    alpha::Any
    x_base::Union{Function,Nothing}
end
# Constructors
function DiscreteDampingSolver(ls::LinearSolver, alpha::Any ; x_base = nothing)
    return DiscreteDampingSolver(ls, alpha, x_base)
end
# Function Definitions
struct DiscreteDampingSolverSymbolicSetup <: SymbolicSetup
    ls::DiscreteDampingSolver
    ss::SymbolicSetup
end
mutable struct DiscreteDampingSolverNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
    ss::DiscreteDampingSolverSymbolicSetup
    A::T
    ns::NumericalSetup
    cache
    x_new::AbstractVector
end
function get_solver_caches(A,solver::DiscreteDampingSolver)
    α_arr = allocate_in_domain(A);
    α_arr .= get_free_dof_values(solver.alpha)
    x_base_arr = isnothing(solver.x_base) ? nothing : fill!(allocate_in_domain(A), 0)
    return α_arr, x_base_arr
end
function update_solver_caches!(ns::DiscreteDampingSolverNumericalSetup,solver::DiscreteDampingSolver;StageOperator = nothing)
    α_arr, x_base_arr = ns.cache
    if(isnothing(StageOperator))
        if(!isnothing(solver.x_base))  x_base_arr .= get_free_dof_values(solver.x_base) end
    else
        if(!isnothing(solver.x_base))  x_base_arr .= get_free_dof_values(solver.x_base(StageOperator.tx)) end
    end
    ns.cache = α_arr, x_base_arr
end
function Gridap.Algebra.symbolic_setup(ls::DiscreteDampingSolver, mat::AbstractMatrix)
    DiscreteDampingSolverSymbolicSetup(ls, Gridap.Algebra.symbolic_setup(ls.ls, mat))
end
function Gridap.Algebra.numerical_setup(
    ss::DiscreteDampingSolverSymbolicSetup,
    A::AbstractMatrix
)    
    x_new = allocate_in_domain(A); # Block Solve Allocation
    fill!(x_new, 0)
    cache = get_solver_caches(A,ss.ls) # Introducing cache
    DiscreteDampingSolverNumericalSetup(ss, A, Gridap.Algebra.numerical_setup(ss.ss, A),cache,x_new)
end
function Gridap.Algebra.numerical_setup!(
    ns::DiscreteDampingSolverNumericalSetup,
    A::AbstractMatrix,
)   
    fill!(ns.x_new, 0.0) # reinitializing solution vector
    ns.A = A
    Gridap.Algebra.numerical_setup!(ns.ns, A)
end
function Gridap.Algebra.solve!(
    x::AbstractVector,
    ns::DiscreteDampingSolverNumericalSetup,
    b::AbstractVector;
    lop::Any = nothing
)   
    solve!(ns.x_new, ns.ns, b)
    copy!(x, ns.x_new)
    update_solver_caches!(ns,ns.ss.ls;StageOperator = lop) # Updating cache
    α_arr, x_base_arr = ns.cache
    if isnothing(ns.ss.ls.x_base)
        x .= α_arr.*x # Damping with zero x_base
    elseif isnothing(lop)
        x .= α_arr.*x + (1 .- α_arr).*x_base_arr # Non transient damping with non zero x_base
    elseif isa(lop, Gridap.ODEs.LinearStageOperator)
        x .= α_arr.*x + (1 .- α_arr).*(x_base_arr .- lop.usx[1])./lop.ws[1] # Transient damping with non zero x_base (Linear Stage Operator)
    elseif isa(lop, Gridap.ODEs.NonLinearStageOperator)
        x .= α_arr.*x + (1 .- α_arr).*(x_base_arr .- lop.usx[1](x))./lop.ws[1] # Transient damping with non zero x_base (Non Linear Stage Operator)
    else
        error("Stage Operator type not recognized for transient damping")
    end
end
# Solver for GridapODE
function Gridap.Algebra.solve!(
    x::AbstractVector,
    ls::DiscreteDampingSolver,
    lop::Gridap.ODEs.LinearStageOperator,
    ns::Nothing,
)
    J = lop.J
    ss = Gridap.Algebra.symbolic_setup(ls, J)
    ns = Gridap.Algebra.numerical_setup(ss, J)
    r = lop.r
    rmul!(r, -1)
    solve!(x, ns, r; lop = lop) # Passing lop for transient damping
    ns
end
function Gridap.Algebra.solve!(
    x::AbstractVector,
    ls::DiscreteDampingSolver,
    lop::Gridap.ODEs.LinearStageOperator,
    ns,
)
    if !lop.reuse
        J = lop.J
        Gridap.Algebra.numerical_setup!(ns, J)
    end
    r = lop.r
    rmul!(r, -1)
    solve!(x, ns, r; lop = lop) # Passing lop for transient damping
    ns
end
# End damping solver setup
##########################

##########################
##### Triangulations #####
##########################

function get_triangulations(model,mem_tag)
    labels = get_face_labeling(model)
    cells_from_nodes_from_solid_faces = Vector{Int32}[]; cells_from_nodes_from_fs_faces_wo_solid = Vector{Int32}[];
    map(local_views(model), local_views(labels)) do model, labels
        # Get the grid topology 
        topo = Gridap.Geometry.get_grid_topology(model)
        # Get the dimension of the model
        D = Gridap.Geometry.num_cell_dims(topo)
        # Extracting all the cells from the "Membrane tag"
        s_Γmask = get_face_mask(labels,mem_tag,D-1)
        s_indices = findall(s_Γmask)
        nodes_from_solid_faces = unique(vcat(get_faces(topo, D-1, 0)[s_indices]...))
        cells_from_nodes_from_solid_faces = unique(vcat(get_faces(topo, 0, D)[nodes_from_solid_faces]...))
    end
    # Creating the triangulations
    ΩΓ_s = GridapDistributed.Geometry.Triangulation(model, Vector{Int32}(cells_from_nodes_from_solid_faces))
    return ΩΓ_s
end

function get_triangulations(model,fs_tag,mem_tag)
    labels = get_face_labeling(model)
    cells_from_nodes_from_solid_faces_wo_fs = Vector{Int32}[]; cells_from_nodes_from_fs_faces_wo_solid = Vector{Int32}[];
    map(local_views(model), local_views(labels)) do model, labels
        # Get the grid topology 
        topo = Gridap.Geometry.get_grid_topology(model)
        # Get the dimension of the model
        D = Gridap.Geometry.num_cell_dims(topo)
        # Extracting all the cells from the "Membrane tag"
        s_Γmask = get_face_mask(labels,mem_tag,D-1)
        s_indices = findall(s_Γmask)
        nodes_from_solid_faces = unique(vcat(get_faces(topo, D-1, 0)[s_indices]...))
        cells_from_nodes_from_solid_faces = unique(vcat(get_faces(topo, 0, D)[nodes_from_solid_faces]...))
        # Extracting all the cells from the "Freesurface tag" and not in the "Membrane tag"
        fs_mask = get_face_mask(labels,fs_tag,D-1)
        fs_indices = findall(fs_mask)
        nodes_from_fs_faces = unique(vcat(get_faces(topo, D-1, 0)[fs_indices]...))
        nodes_from_fs_faces_wo_solid = setdiff(nodes_from_fs_faces, nodes_from_solid_faces)
        cells_from_nodes_from_fs_faces_wo_solid = unique(vcat(get_faces(topo, 0, D)[nodes_from_fs_faces_wo_solid]...))
        # Extracting all the cells from the "Membrane tag" and not in the "Freesurface tag"
        cells_from_nodes_from_solid_faces_wo_fs = setdiff(cells_from_nodes_from_solid_faces, cells_from_nodes_from_fs_faces_wo_solid)
    end
    # Creating the triangulations
    ΩΓ_fs = GridapDistributed.Geometry.Triangulation(model, Vector{Int32}(cells_from_nodes_from_fs_faces_wo_solid))
    ΩΓ_s = GridapDistributed.Geometry.Triangulation(model, Vector{Int32}(cells_from_nodes_from_solid_faces_wo_fs))
    return ΩΓ_fs, ΩΓ_s
end

function select_triangulation(a::GridapDistributed.DistributedTriangulation,b::GridapDistributed.DistributedTriangulation)
  @assert a.model === b.model
  trians = map(a.trians,b.trians) do a,b
    InterfaceTriangulation(a,b).plus
  end
  GridapDistributed.DistributedTriangulation(trians,a.model)
end

function get_interface_triangulation_side(a, b , side)
    @assert a.model === b.model
    trians = map(a.trians,b.trians) do a,b
        getfield(InterfaceTriangulation(a,b),side)
    end
    GridapDistributed.DistributedTriangulation(trians,a.model)
end

# Macro to conditionally execute code for vtk output
macro vtk(expr)
    return esc(quote
        if vtkoutput
            $expr
        end
    end)
end

# Notation for scientific string formatting
function sci_str(x::Real)
    x == 0 && return "0.00e0"
    exponent = floor(Int, log10(abs(x)))
    mantissa = round(x / 10^exponent, digits=2)
    return "$(mantissa)e$(exponent)"
end

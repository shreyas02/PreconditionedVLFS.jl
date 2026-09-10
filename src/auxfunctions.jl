################################
# Begin DiscreteDampingSolver setup
struct DiscreteDampingSolver <: LinearSolver
  ls::LinearSolver
  alpha::Any
  x_base::Union{Function,Nothing}
end

# Constructors
function DiscreteDampingSolver(ls::LinearSolver, alpha::Any; x_base = nothing)
  return DiscreteDampingSolver(ls, alpha, x_base)
end

# Function Definitions
struct DiscreteDampingSolverSymbolicSetup <: SymbolicSetup
  ls::DiscreteDampingSolver
  ss::SymbolicSetup
end

mutable struct DiscreteDampingSolverNumericalSetup{
  T<:AbstractMatrix,
} <: NumericalSetup
  ss::DiscreteDampingSolverSymbolicSetup
  A::T
  ns::NumericalSetup
  cache
  x_new::AbstractVector
end

function get_solver_caches(A, solver::DiscreteDampingSolver)
  α_arr = allocate_in_domain(A)
  α_arr .= get_free_dof_values(solver.alpha)
  x_base_arr = if isnothing(solver.x_base)
    nothing
  else
    fill!(allocate_in_domain(A), 0)
  end
  return α_arr, x_base_arr
end

function update_solver_caches!(
  ns::DiscreteDampingSolverNumericalSetup,
  solver::DiscreteDampingSolver;
  stage_operator = nothing,
)
  α_arr, x_base_arr = ns.cache
  if isnothing(stage_operator)
    if !isnothing(solver.x_base)
      x_base_arr .= get_free_dof_values(solver.x_base)
    end
  else
    if !isnothing(solver.x_base)
      x_base_arr .= get_free_dof_values(solver.x_base(stage_operator.tx))
    end
  end
  ns.cache = α_arr, x_base_arr
end

function Gridap.Algebra.symbolic_setup(
  ls::DiscreteDampingSolver,
  mat::AbstractMatrix,
)
  return DiscreteDampingSolverSymbolicSetup(
    ls,
    Gridap.Algebra.symbolic_setup(ls.ls, mat),
  )
end

function Gridap.Algebra.numerical_setup(
  ss::DiscreteDampingSolverSymbolicSetup,
  A::AbstractMatrix,
)
  x_new = allocate_in_domain(A) # Block Solve Allocation
  fill!(x_new, 0)
  cache = get_solver_caches(A, ss.ls) # Introducing cache
  return DiscreteDampingSolverNumericalSetup(
    ss,
    A,
    Gridap.Algebra.numerical_setup(ss.ss, A),
    cache,
    x_new,
  )
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
  lop::Any = nothing,
)
  solve!(ns.x_new, ns.ns, b)
  copy!(x, ns.x_new)
  update_solver_caches!(ns, ns.ss.ls; stage_operator = lop) # Updating cache
  α_arr, x_base_arr = ns.cache
  if isnothing(ns.ss.ls.x_base)
    x .= α_arr .* x # Damping with zero x_base
  elseif isnothing(lop)
    x .= α_arr .* x + (1 .- α_arr) .* x_base_arr
  elseif isa(lop, Gridap.ODEs.LinearStageOperator)
    x .= (
      α_arr .* x +
      (1 .- α_arr) .* (x_base_arr .- lop.usx[1]) ./ lop.ws[1]
    )
  elseif isa(lop, Gridap.ODEs.NonLinearStageOperator)
    x .= (
      α_arr .* x +
      (1 .- α_arr) .* (x_base_arr .- lop.usx[1](x)) ./ lop.ws[1]
    )
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
  return ns
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
  return ns
end

# End damping solver setup
##########################

##########################
##### Triangulations #####
##########################

function get_triangulations(model, mem_tag)
  labels = get_face_labeling(model)
  cells_from_nodes_from_solid_faces = Vector{Int32}[]
  map(local_views(model), local_views(labels)) do model, labels
    # Get the grid topology
    topo = Gridap.Geometry.get_grid_topology(model)
    # Get the dimension of the model
    D = Gridap.Geometry.num_cell_dims(topo)
    # Extracting all the cells from the "Membrane tag"
    s_Γmask = get_face_mask(labels, mem_tag, D - 1)
    s_indices = findall(s_Γmask)
    nodes_from_solid_faces = unique(
      vcat(get_faces(topo, D - 1, 0)[s_indices]...),
    )
    cells_from_nodes_from_solid_faces = unique(
      vcat(get_faces(topo, 0, D)[nodes_from_solid_faces]...),
    )
  end
  # Creating the triangulations
  ΩΓ_s = GridapDistributed.Geometry.Triangulation(
    model,
    Vector{Int32}(cells_from_nodes_from_solid_faces),
  )
  return ΩΓ_s
end

function get_triangulations(model, fs_tag, mem_tag)
  labels = get_face_labeling(model)
  cells_from_nodes_from_solid_faces_wo_fs = Vector{Int32}[]
  cells_from_nodes_from_fs_faces_wo_solid = Vector{Int32}[]
  map(local_views(model), local_views(labels)) do model, labels
    # Get the grid topology
    topo = Gridap.Geometry.get_grid_topology(model)
    # Get the dimension of the model
    D = Gridap.Geometry.num_cell_dims(topo)
    # Extracting all the cells from the "Membrane tag"
    s_Γmask = get_face_mask(labels, mem_tag, D - 1)
    s_indices = findall(s_Γmask)
    nodes_from_solid_faces = unique(
      vcat(get_faces(topo, D - 1, 0)[s_indices]...),
    )
    cells_from_nodes_from_solid_faces = unique(
      vcat(get_faces(topo, 0, D)[nodes_from_solid_faces]...),
    )
    # Extracting cells from the "Freesurface tag" and not the "Membrane tag"
    fs_mask = get_face_mask(labels, fs_tag, D - 1)
    fs_indices = findall(fs_mask)
    nodes_from_fs_faces = unique(
      vcat(get_faces(topo, D - 1, 0)[fs_indices]...),
    )
    nodes_from_fs_faces_wo_solid = setdiff(
      nodes_from_fs_faces,
      nodes_from_solid_faces,
    )
    cells_from_nodes_from_fs_faces_wo_solid = unique(
      vcat(get_faces(topo, 0, D)[nodes_from_fs_faces_wo_solid]...),
    )
    # Extracting cells from the "Membrane tag" and not the "Freesurface tag"
    cells_from_nodes_from_solid_faces_wo_fs = setdiff(
      cells_from_nodes_from_solid_faces,
      cells_from_nodes_from_fs_faces_wo_solid,
    )
  end
  # Creating the triangulations
  ΩΓ_fs = GridapDistributed.Geometry.Triangulation(
    model,
    Vector{Int32}(cells_from_nodes_from_fs_faces_wo_solid),
  )
  ΩΓ_s = GridapDistributed.Geometry.Triangulation(
    model,
    Vector{Int32}(cells_from_nodes_from_solid_faces_wo_fs),
  )
  return ΩΓ_fs, ΩΓ_s
end

function select_triangulation(
  a::GridapDistributed.DistributedTriangulation,
  b::GridapDistributed.DistributedTriangulation,
)
  @assert a.model === b.model
  trians = map(a.trians, b.trians) do a, b
    InterfaceTriangulation(a, b).plus
  end
  return GridapDistributed.DistributedTriangulation(trians, a.model)
end

function get_interface_triangulation_side(a, b, side)
  @assert a.model === b.model
  trians = map(a.trians, b.trians) do a, b
    getfield(InterfaceTriangulation(a, b), side)
  end
  return GridapDistributed.DistributedTriangulation(trians, a.model)
end

# Macro to conditionally execute code for vtk output
macro vtk(expr)
  return esc(
    quote
      if vtkoutput
        $expr
      end
    end,
  )
end

# Notation for scientific string formatting
function sci_str(x::Real)
  x == 0 && return "0.00e0"
  exponent = floor(Int, log10(abs(x)))
  mantissa = round(x / 10^exponent, digits = 2)
  return "$(mantissa)e$(exponent)"
end

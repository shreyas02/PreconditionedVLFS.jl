# SolverWrapper setup
struct SolverWrapper <: LinearSolver
  ls::LinearSolver
end

struct SolverWrapperSymbolicSetup <: SymbolicSetup
  ss::SymbolicSetup
end

struct SolverWrapperNumericalSetup <: NumericalSetup
  ns::NumericalSetup
  x_new::AbstractVector
end

function Gridap.Algebra.symbolic_setup(ls::SolverWrapper, A::AbstractMatrix)
  return SolverWrapperSymbolicSetup(symbolic_setup(ls.ls, A))
end

function Gridap.Algebra.numerical_setup(ss::SolverWrapperSymbolicSetup, A::AbstractMatrix)
  x_new = allocate_in_domain(A)
  fill!(x_new, 0)
  return SolverWrapperNumericalSetup(numerical_setup(ss.ss, A), x_new)
end

function Gridap.Algebra.numerical_setup!(ns::SolverWrapperNumericalSetup, A::AbstractMatrix)
  numerical_setup!(ns.ns, A)
  return ns
end

function Gridap.Algebra.solve!(x::AbstractVector, ns::SolverWrapperNumericalSetup, b::AbstractVector)
  fill!(ns.x_new, 0)
  solve!(ns.x_new, ns.ns, b)
  copy!(x, ns.x_new)
  return x
end

"""
    safeintegral(f, dΣ)

Integrate a distributed surface field, skipping empty local integration domains.
This avoids evaluating incompatible placeholder bases on empty surface partitions.
"""
function safeintegral(f, dΣ)
  contributions = map(local_views(f), local_views(dΣ)) do f_loc, measure_loc
    quad = Gridap.CellData.get_cell_quadrature(measure_loc)
    if num_cells(get_triangulation(quad)) == 0
      return Gridap.CellData.DomainContribution()
    end
    return ∫(f_loc) * measure_loc
  end
  return GridapDistributed.DistributedDomainContribution(contributions)
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

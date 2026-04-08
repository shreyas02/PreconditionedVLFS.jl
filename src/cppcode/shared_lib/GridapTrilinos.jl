# Load the module and generate the functions
module GridapTrilinos
  using CxxWrap
  @wrapmodule(() -> joinpath(@__DIR__, "libmpicpp"))

  # cpp structs and functions to be exported to julia
  export TrilinosParallel, KokkosInitialize, KokkosFinalize
  export SolverResult, SolverResultAllocated, SolverResultDereferenced
  # properties of the SolverResult struct to be exported to julia
  export num_iters, residual, solve_time, name, residuals, verbose, depth

  # cpp vector to julia vector 
  residuals(result::SolverResult) = collect(_residualsCxx(result))

  const _solver_result_properties = (:num_iters, :residual, :solve_time, :residuals, :name, :verbose, :depth)

  function Base.getproperty(result::Union{SolverResultAllocated,SolverResultDereferenced}, element::Symbol)
    if element === :num_iters
      return () -> num_iters(result)
    elseif element === :residual
      return () -> residual(result)
    elseif element === :solve_time
      return () -> solve_time(result)
    elseif element === :residuals
      return () -> residuals(result)
    elseif element === :name
      return () -> name(result)
    elseif element === :verbose
      return () -> verbose(result)
    elseif element === :depth
      return () -> depth(result)
    end
    return getfield(result, element)
  end

  Base.propertynames(::Union{SolverResultAllocated,SolverResultDereferenced}, private::Bool=false) =
    private ? (_solver_result_properties..., fieldnames(SolverResultAllocated)...) : _solver_result_properties

  function __init__()
    @initcxx
    KokkosInitialize()
    atexit(KokkosFinalize)
  end
end

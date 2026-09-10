module PreconditionedVLFS

# Header code
include("julia_headers.jl")

# Exporting the package elements

# Problems and their parameters
export wsi2d, WSI2D_params
export wsi3d, WSI3D_params
export toyrichardson, ToyRichardsonParams
export periodic2D, Periodic2D_params
export periodic3D, Periodic3D_params

# Auxiliary functions
export TrilinosSolve
export DiscreteDampingSolver
export get_triangulations
export select_triangulation
export get_interface_triangulation_side
export sci_str

# Gridap Trilinos wrapper
export SolverResult, SolverResultAllocated, SolverResultDereferenced
export KokkosInitialize, KokkosFinalize
export num_iters, residual, solve_time, name, verbose, depth

# Exporting all the macro helpers
export @vtk

# Include auxiliary functions
include("auxfunctions.jl")

# Application source files
include("WSI2D.jl")
include("WSI3D.jl")
include("ToyRichardson.jl")
include("Periodic2D.jl")
include("Periodic3D.jl")

end

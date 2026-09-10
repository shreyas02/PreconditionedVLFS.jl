module Periodic3DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson

# Case map:
#   base_case:      periodic3d_base_case
#   strong_scaling: periodic3d_strong_scaling
#   weak_scaling:   periodic3d_weak_scaling
#
# Each case is stored through DrWatson's `produce_or_load`.  The
# `case` field is included in the generated data filename.

function base_case()
  with_mpi() do distribute
    parts = (8, 2, 1)

    # Function to run the code
    function run_src(case::Periodic3D_params)
      solver_stats = periodic3D(distribute, parts, case)
      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(case)
      )
      return merge(
        config_dict,
        Dict(
          "fluid" => solver_stats.fluid,
          "solid" => solver_stats.solid,
          "outer" => solver_stats.outer,
        ),
      )
    end
    # Base case | Construction of the required parameters
    case = Periodic3D_params(
      nprocs = MPI.Comm_size(MPI.COMM_WORLD),
      rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
      case = "periodic3d_base_case",
      H = 10, # Height of the domain
      Lf = 16 * π, # Length of the entire domain
      Ly = 10, # Width of the domain
      dt = 0.1,
      tF = 0.2,
    )

    # Path for this test case
    path = mkpath("$(datadir("periodic3d", "base_case"))")
    # Run the code
    data, _ = produce_or_load(run_src, case, path)
  end
end

function strong_scaling()
  with_mpi() do distribute
    parts = (MPI.Comm_size(MPI.COMM_WORLD) ÷ 4, 4, 1)

    # Function to run the code
    function run_src(case::Periodic3D_params)
      solver_stats = periodic3D(distribute, parts, case)
      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(case)
      )
      return merge(
        config_dict,
        Dict(
          "fluid" => solver_stats.fluid,
          "solid" => solver_stats.solid,
          "outer" => solver_stats.outer,
        ),
      )
    end
    # Strong scaling test case | Construction of the required parameters
    # We run a single case for 10 iterations and store the solver timings.
    for i = 1:10
      case = Periodic3D_params(
        nprocs = MPI.Comm_size(MPI.COMM_WORLD),
        rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
        iter = i,
        case = "periodic3d_strong_scaling",
        H = 10, # Height of the domain
        # Length for 64 processes, kept fixed for strong scaling.
        Lf = 64 * π,
        Ly = 10, # Width of the domain
        dt = 0.1,
        tF = 0.2,
      )

      # Path for this test case
      path = mkpath("$(datadir("periodic3d", "strong_scaling"))")
      # Run the code
      data, _ = produce_or_load(run_src, case, path)
    end
  end
end

function weak_scaling()
  with_mpi() do distribute
    parts = (MPI.Comm_size(MPI.COMM_WORLD) ÷ 2, 2, 1)

    # Function to run the code
    function run_src(case::Periodic3D_params)
      solver_stats = periodic3D(distribute, parts, case)
      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(case)
      )
      return merge(
        config_dict,
        Dict(
          "fluid" => solver_stats.fluid,
          "solid" => solver_stats.solid,
          "outer" => solver_stats.outer,
        ),
      )
    end
    # Weak scaling test case | Construction of the required parameters
    # We run a single case for 10 iterations and store the solver timings.
    for i = 1:10
      case = Periodic3D_params(
        nprocs = MPI.Comm_size(MPI.COMM_WORLD),
        rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
        case = "periodic3d_weak_scaling",
        iter = i,
        H = 10, # Height of the domain
        Lf = 0.25 * π * MPI.Comm_size(MPI.COMM_WORLD),
        Ly = 10, # Width of the domain
        dt = 0.1,
        tF = 0.2,
      )

      # Path for this test case
      path = mkpath("$(datadir("periodic3d", "weak_scaling"))")
      # Run the code
      data, _ = produce_or_load(run_src, case, path)
    end
  end
end

end # module

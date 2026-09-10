module Periodic2DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson

# Case map:
#   base_case:      periodic2d_base_case
#   strong_scaling: periodic2d_strong_scaling
#   weak_scaling:   periodic2d_weak_scaling
#
# Each case is stored through DrWatson's `produce_or_load`.  The
# `case` field is included in the generated data filename.

function base_case()
  with_mpi() do distribute
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)

    # Function to run the code
    function run_src(case::Periodic2D_params)
      solver_stats = periodic2D(distribute, parts, case)
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
    case = Periodic2D_params(
      nprocs = MPI.Comm_size(MPI.COMM_WORLD),
      rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
      case = "periodic2d_base_case",
      H = 10, # Height of the domain
      Lf = 20 * π, # Length of the entire domain
      hs = 0.01, # Floating structure thickness
      dt = 0.1,
      tF = 100.0,
    )

    # Path for this test case
    path = mkpath("$(datadir("periodic2d", "base_case"))")
    # Run the code
    data, _ = produce_or_load(run_src, case, path)
  end
end

function strong_scaling()
  with_mpi() do distribute
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)

    # Function to run the code
    function run_src(case::Periodic2D_params)
      solver_stats = periodic2D(distribute, parts, case)
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
      case = Periodic2D_params(
        nprocs = MPI.Comm_size(MPI.COMM_WORLD),
        rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
        iter = i,
        case = "periodic2d_strong_scaling",
        H = 10, # Height of the domain
        # Length of the entire domain for the smallest case with 64 ranks.
        # Used for all cases.
        Lf = 360 * π * 2,
        hs = 0.01, # Floating structure thickness
        dt = 0.1,
        tF = 1.2,
      )

      # Path for this test case
      path = mkpath("$(datadir("periodic2d", "strong_scaling"))")
      # Run the code
      data, _ = produce_or_load(run_src, case, path)
    end
  end
end

function weak_scaling()
  with_mpi() do distribute
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)

    # Function to run the code
    function run_src(case::Periodic2D_params)
      solver_stats = periodic2D(distribute, parts, case)
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
      case = Periodic2D_params(
        nprocs = MPI.Comm_size(MPI.COMM_WORLD),
        rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
        case = "periodic2d_weak_scaling",
        iter = i,
        H = 10, # Height of the domain
        # Length of the entire domain for the smallest case with 64 ranks.
        Lf = 11.25 * π * MPI.Comm_size(MPI.COMM_WORLD),
        hs = 0.01, # Floating structure thickness
        dt = 0.1,
        tF = 1.2,
      )

      # Path for this test case
      path = mkpath("$(datadir("periodic2d", "weak_scaling"))")
      # Run the code
      data, _ = produce_or_load(run_src, case, path)
    end
  end
end

end # module

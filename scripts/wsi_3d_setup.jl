module WSI3DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson, TimerOutputs
using Plots, DataFrames, Roots

include("mesh_wsi_3d_1.jl")
using .WSI3DMesh1

# A simple case that can be run easily on a laptop.

function case_1()
  with_mpi() do distribute

    # Generate ranks
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)
    ranks = distribute(LinearIndices((prod(parts),)))

    # Function to run the code
    function run_src(params::WSI3D_params)
      solver_stats = wsi3d(distribute, parts, params)

      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(params)
      )

      return merge(
        config_dict,
        Dict(
          "fluid" => solver_stats.fluid,
          "solid" => solver_stats.solid,
          "freesurface" => solver_stats.freesurface,
          "outer" => solver_stats.outer,
        ),
      )
    end

    # Construction of the required parameters

    # Case number
    case_name = "case_1"

    # Geometric parameters
    H = 1.1
    Lm = pi
    Lf = 9 * pi
    Ly = pi
    hs = 0.01

    # Damping parameters
    Lfd = 3 * pi
    Lfd1 = 0.5 * pi
    Ld = 7.5 * pi
    Ld1 = Lf - 0.5

    # Temporal parameters
    ρ∞ = 0.5
    t0 = 0.0
    tF = 20.0
    dt = 0.1

    # Physical parameters
    ρf = 1000.0 # Fluid density
    ρs = 100 # Solid density
    g = 9.81 # Acceleration due to gravity
    T = 0.9 * ρf * g # Solid stiffness parameter

    # Wave parameters
    kλ = 3.0 # Wave number
    ω = sqrt(g * kλ * tanh(kλ * H)) # Wave frequency in radians
    η₀ = 0.01 # surface elevation
    ϕ = 0 # wave phase difference

    # Mesh generations
    meshpath = WSI3DMesh1.create_mesh(ranks, H)

    # Post-processing parameters
    vtkoutput = true

    case = WSI3D_params(
      # MPI parameters and case name
      nprocs = MPI.Comm_size(MPI.COMM_WORLD),
      rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
      case = case_name,

      # Geometric parameters
      H = H,
      Lm = Lm,
      Ly = Ly,
      Lf = Lf,
      hs = hs,
      meshpath = meshpath,

      # Damping parameters
      Lfd = Lfd,
      Lfd1 = Lfd1,
      Ld = Ld,
      Ld1 = Ld1,

      # Temporal parameters
      ρ∞ = ρ∞,
      t0 = t0,
      tF = tF,
      dt = dt,

      # Physical parameters
      ρf = ρf, # Fluid density
      ρs = ρs, # Solid density
      g = g, # Acceleration due to gravity
      T = T, # Solid stiffness parameter

      # Wave parameters
      kλ = kλ,
      η₀ = η₀,
      ω = ω,
      ϕ = ϕ,

      # Post-processing parameters
      vtkoutput = vtkoutput,
    )

    path = mkpath("$(datadir("wsi_3d", "case_1"))")
    filename = savename(case; ignores = [:meshpath])

    produce_or_load(run_src, case, path; filename = filename)
  end
end

end

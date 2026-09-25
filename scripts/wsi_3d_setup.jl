module WSI3DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson, TimerOutputs

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

    # Non dimensionalization parameters
    Lref = pi # Physical Membrane length
    g = 9.81 # Acceleration due to gravity
    Uref = sqrt(g * Lref) # Reference velocity
    Tref = Lref / Uref # Reference time

    # Geometric parameters
    H = 1.1/Lref
    Lm = 1
    Lf = 9 * pi/Lref
    Ly = pi/Lref
    hs = 0.01/Lref

    # Damping parameters
    Lfd = 3 * pi/Lref
    Lfd1 = 0.5 * pi/Lref
    Ld = 7.5 * pi/Lref
    Ld1 = Lf - 0.5/Lref

    # Temporal parameters
    ρ∞ = 0.5
    t0 = 0.0/Tref
    tF = 20.0/Tref
    dt = 0.1/Tref

    # Physical parameters.
    M = 0.045 # Non dimensionalized reduced mass parameter
    τ = 0.025 # Non dimensional pretension parameter

    # Wave parameters
    kλ_dim = 3.0;  kλ = kλ_dim * Lref # Wave number
    ω_dim = sqrt(g * kλ_dim * tanh(kλ_dim * H * Lref)); ω = ω_dim * Tref # Wave frequency in radians
    η₀ = 0.01/Lref # surface elevation
    ϕ = 0 # wave phase difference

    # Mesh generations
    meshpath = WSI3DMesh1.create_mesh(ranks, Lref)

    # Post-processing parameters
    vtkoutput = true

    case = WSI3D_params(
      # MPI parameters and case name
      nprocs = MPI.Comm_size(MPI.COMM_WORLD),
      rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
      case = case_name,

      # Reference dimensional state
      Lref = Lref,
      Tref = Tref,

      # Physical parameters
      M = M,
      τ = τ,

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

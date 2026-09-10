module WSI2DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson
using TimerOutputs
using Roots

include("mesh_wsi_2d_1.jl")
using .WSI2DMesh1

# 2D wave-structure interaction validation cases based on "Agarwal, Shagun,
# et al. “Dynamic Analysis of Viscoelastic Floating Membranes Using Monolithic
# Finite Element Method.” Journal of Fluids and Structures, vol. 129, Oct.
# 2024, p. 104167, https://doi.org/10.1016/j.jfluidstructs.2024.104167.
# Accessed 10 Sept. 2026."

function case_1()
  with_mpi() do distribute
    # One-dimensional MPI decomposition used by the 2D mesh generator and
    # solver.
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)
    ranks = distribute(LinearIndices((prod(parts),)))

    # DrWatson callback: run the solver and persist the full parameter set
    # together with the nonlinear outer-iteration history.
    function run_src(params::WSI2D_params)
      solver_stats = wsi2d(distribute, parts, params)
      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(params)
      )

      return merge(
        config_dict,
        Dict(
          "outer" => solver_stats.outer,
        ),
      )
    end

    case_name = "case_1"

    # Geometry: membrane centered in the physical wave tank, with damping
    # regions added upstream and downstream.
    H = 10.0 # Height of the domain
    domain = 18 * H # Length of the required domain
    Lm = 2 * H # Membrane length
    damp = 15 * H # Inlet and outlet damping length
    Lf = domain + 2 * damp # Total domain length
    hs = 0.01 # Thickness of the membrane
    meshpath = WSI2DMesh1.create_mesh(ranks, H, damp, domain, Lm)

    # Damping-zone extents used by the inlet/outlet absorbing layers.
    Lfd = damp
    Lfd1 = 0.5 * Lm
    Ld = Lf - damp
    Ld1 = Lf - 0.5 * Lm

    # Time integration parameters.
    ρ∞ = 1.0
    t0 = 0.0
    tF = 120.0
    dt = 0.05

    # Physical parameters.
    ρf = 1025.0 # Fluid density
    ρs = (ρf * H * 0.090) / hs # Solid density
    g = 9.81 # Acceleration due to gravity
    T = 0.1 * ρf * g * H * H # Solid stiffness parameter

    # Incident-wave parameters. kλ is computed from the finite-depth
    # dispersion relation for the chosen frequency.
    η₀ = 0.1
    ω = 2.4 # Wave frequency
    kλ = find_zero(kλ -> g * kλ * tanh(kλ * H) - ω^2, (0.01, 10.0))
    ϕ = 0.0

    # Post-processing controls.
    vtkoutput = true

    params = WSI2D_params(
      # MPI parameters and case name
      nprocs = MPI.Comm_size(MPI.COMM_WORLD),
      rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
      case = case_name,

      # Geometric parameters
      H = H,
      Lm = Lm,
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
      ρf = ρf,
      ρs = ρs,
      g = g,
      T = T,

      # Wave parameters
      kλ = kλ,
      η₀ = η₀,
      ω = ω,
      ϕ = ϕ,

      # Post-processing parameters
      vtkoutput = vtkoutput,
    )

    path = mkpath("$(datadir("wsi_2d", "case_1"))")
    filename = savename(params; ignores = [:meshpath])

    # Reuse existing DrWatson output when this exact configuration was already
    # run.
    produce_or_load(run_src, params, path; filename = filename)
  end
end

# Density sensitivity analysis for the case_1 setup.

function case_1_density_sweep()
  with_mpi() do distribute
    # One-dimensional MPI decomposition used by the 2D mesh generator and
    # solver.
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)
    ranks = distribute(LinearIndices((prod(parts),)))

    # DrWatson callback: run the solver and persist the full parameter set
    # together with the nonlinear outer-iteration history.
    function run_src(params::WSI2D_params)
      solver_stats = wsi2d(distribute, parts, params)
      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(params)
      )

      return merge(
        config_dict,
        Dict(
          "outer" => solver_stats.outer,
        ),
      )
    end

    case_name = "case_1_density_sweep"

    # Geometry is fixed across the density sweep so only the density ratio
    # changes.
    H = 10.0 # Height of the domain
    domain = 2.1 * H # Length of the required domain
    Lm = 2 * H # Membrane length
    damp = 5 * H # Inlet and outlet damping length
    Lf = domain + 2 * damp # Total domain length
    hs = 0.01 # Thickness of the membrane
    meshpath = WSI2DMesh1.create_mesh(ranks, H, damp, domain, Lm)

    # Damping-zone extents used by the inlet/outlet absorbing layers.
    Lfd = damp
    Lfd1 = 0.5 * Lm
    Ld = Lf - damp
    Ld1 = Lf - 0.5 * Lm

    # Time integration parameters.
    ρ∞ = 1.0
    t0 = 0.0
    tF = 10.0
    dt = 0.1

    # Physical parameters. The sweep scales the reference membrane density.
    ρf = 1025.0 # Fluid density
    ρs_standard = (ρf * H * 0.090) / hs # Solid density
    density_ratios = [0.0001, 0.001, 0.01, 0.1]
    ρs_sweep = density_ratios .* ρs_standard # Sweep over solid density values
    g = 9.81 # Acceleration due to gravity
    T = 0.9 * ρf * g # * ρf * g * H * H # Solid stiffness parameter

    # Incident-wave parameters. kλ is computed from the finite-depth
    # dispersion relation for the chosen frequency.
    η₀ = 0.1
    ω = 2.0 # Wave frequency
    kλ = find_zero(kλ -> g * kλ * tanh(kλ * H) - ω^2, (0.01, 10.0))
    ϕ = 0.0

    # Post-processing controls and in-memory sweep summary.
    vtkoutput = true
    density_sweep_data = []

    for (density_ratio, ρs) in zip(density_ratios, ρs_sweep)
      params = WSI2D_params(
        # MPI parameters and case name
        nprocs = MPI.Comm_size(MPI.COMM_WORLD),
        rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
        case = case_name,

        # Geometric parameters
        H = H,
        Lm = Lm,
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
        ρf = ρf,
        ρs = ρs,
        g = g,
        T = T,

        # Wave parameters
        kλ = kλ,
        η₀ = η₀,
        ω = ω,
        ϕ = ϕ,

        # Post-processing parameters
        vtkoutput = vtkoutput,
      )

      path = mkpath("$(datadir("wsi_2d", case_name))")
      filename = savename(
        Dict(
          :case => case_name,
          :density_ratio => density_ratio,
          :nprocs => params.nprocs,
          :rank => params.rank,
        ),
      )

      data, _ = produce_or_load(run_src, params, path; filename = filename)
      push!(
        density_sweep_data,
        (solid_fluid_density_ratio = ρs / ρf, data = data),
      )
    end # density sweep

    # MPI barrier through a distributed reduction before leaving the case.
    sum(ranks)
  end # MPI context
end # case_1_density_sweep

# Membrane-length sensitivity analysis for the case_1 setup.

function case_1_length_sweep()
  with_mpi() do distribute
    # One-dimensional MPI decomposition used by the 2D mesh generator and
    # solver.
    parts = (MPI.Comm_size(MPI.COMM_WORLD), 1)
    ranks = distribute(LinearIndices((prod(parts),)))

    # DrWatson callback: run the solver and persist the full parameter set
    # together with the nonlinear outer-iteration history.
    function run_src(params::WSI2D_params)
      solver_stats = wsi2d(distribute, parts, params)
      config_dict = Dict(
        string(k) => v for (k, v) in DrWatson.struct2dict(params)
      )

      return merge(
        config_dict,
        Dict(
          "outer" => solver_stats.outer,
        ),
      )
    end

    case_name = "case_1_length_sweep"

    # Base geometry. Each membrane length gets its own mesh inside the loop.
    H = 10.0 # Height of the domain
    length_ratios = [1.0, 5.0, 10.0, 100.0]
    Lm_sweep = length_ratios .* H # Sweep over membrane lengths
    hs = 0.01 # Thickness of the membrane

    # Time integration parameters.
    ρ∞ = 1.0
    t0 = 0.0
    tF = 10.0 # 2 seconds to reach steady state, 7 seconds to record data
    dt = 0.1

    # Physical parameters held fixed across the length sweep.
    ρf = 1025.0 # Fluid density
    ρs = (ρf * H * 0.090) / hs # Solid density
    g = 9.81 # Acceleration due to gravity
    T = 0.9 * ρf * g # * ρf * g * H * H # Solid stiffness parameter

    # Incident-wave parameters. kλ is computed from the finite-depth
    # dispersion relation for the chosen frequency.
    η₀ = 0.1
    ω = 2.0 # Wave frequency
    kλ = find_zero(kλ -> g * kλ * tanh(kλ * H) - ω^2, (0.01, 10.0))
    print("Calculated wave number kλ: ", kλ)
    ϕ = 0.0

    # Post-processing controls and in-memory sweep summary.
    vtkoutput = true
    length_sweep_data = []

    for (length_ratio, Lm) in zip(length_ratios, Lm_sweep)
      damp = 5 * H # Inlet and outlet damping length

      # Keep a small clearance around the membrane while varying its length.
      domain = Lm + 0.2 # Length of the required domain
      Lf = domain + 2 * damp # Total domain length

      # Damping-zone extents for this mesh.
      Lfd = damp
      Ld = Lf - damp
      Lfd1 = 0.5 * Lm
      Ld1 = Lf - 0.5 * Lm

      meshpath = WSI2DMesh1.create_mesh(ranks, H, damp, domain, Lm)

      params = WSI2D_params(
        # MPI parameters and case name
        nprocs = MPI.Comm_size(MPI.COMM_WORLD),
        rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
        case = case_name,

        # Geometric parameters
        H = H,
        Lm = Lm,
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
        ρf = ρf,
        ρs = ρs,
        g = g,
        T = T,

        # Wave parameters
        kλ = kλ,
        η₀ = η₀,
        ω = ω,
        ϕ = ϕ,

        # Post-processing parameters
        vtkoutput = vtkoutput,
      )

      path = mkpath("$(datadir("wsi_2d", case_name))")
      filename = savename(
        Dict(
          :case => case_name,
          :length_ratio => length_ratio,
          :Lm => Lm,
          :nprocs => params.nprocs,
          :rank => params.rank,
        ),
      )

      data, _ = produce_or_load(run_src, params, path; filename = filename)
      push!(length_sweep_data, (membrane_length = Lm, data = data))
    end # length sweep

    # MPI barrier through a distributed reduction before leaving the case.
    sum(ranks)
  end # MPI context
end # case_1_length_sweep

end

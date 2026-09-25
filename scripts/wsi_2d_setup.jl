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

    # Non dimensionalization parameters
    Lref = 20.0 # Physical Membrane length
    g = 9.81 # Acceleration due to gravity
    Uref = sqrt(g * Lref) # Reference velocity
    Tref = Lref / Uref # Reference time

    # Non dimensionalized parameters
    # Geometry: membrane centered in the physical wave tank, with damping
    # regions added upstream and downstream.
    Lm = 1.0 # Membrane length
    H = Lm / 2 # Height of the domain
    domain = 9 * Lm # Length of the required domain
    damp = 7.5 * Lm # Inlet and outlet damping length
    Lf = domain + 2 * damp # Total domain length
    hs = 0.01 / Lref # Thickness of the membrane
    meshpath = WSI2DMesh1.create_mesh(ranks, H, damp, domain, Lm)

    # Damping-zone extents used by the inlet/outlet absorbing layers.
    Lfd = damp
    Lfd1 = 0.5 * Lm
    Ld = Lf - damp
    Ld1 = Lf - 0.5 * Lm

    # Time integration parameters.
    ρ∞ = 0.5
    t0 = 0.0 / Tref
    tF = 120.0 / Tref
    dt = 0.1 / Tref

    # Physical parameters.
    M = 0.045 # Non dimensionalized reduced mass parameter
    τ = 0.025 # Non dimensional pretension parameter

    # Incident-wave parameters
    η₀_dim = 0.1 # Wave amplitude [m]
    ω_dim = 2.0 # Angular frequency [rad/s]
    # Nondimensional wave parameters
    η₀ = η₀_dim / Lref
    ω = ω_dim * Tref
    # Wavenumber from dimensional finite-depth dispersion relation
    kλ_dim = find_zero(
      k -> g * k * tanh(k * H * Lref) - ω_dim^2,
      (0.01, 10.0),
    )
    kλ = kλ_dim * Lref # Nondimensional wavenumber
    ϕ = 0.0

    # Post-processing controls.
    vtkoutput = true

    params = WSI2D_params(
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
      Lf = Lf,
      Lm = Lm,
      hs = hs,
      meshpath = meshpath,

      # Damping parameters
      Lfd = Lfd,
      Lfd1 = Lfd1,
      Ld = Ld,
      Ld1 = Ld1,

      # Wave parameters
      kλ = kλ,
      η₀ = η₀,
      ϕ = ϕ,
      ω = ω,

      # Temporal parameters
      ρ∞ = ρ∞,
      dt = dt,
      t0 = t0,
      tF = tF,

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

    # Non dimensionalization parameters
    Lref = 20.0 # Physical Membrane length
    g = 9.81 # Acceleration due to gravity
    Uref = sqrt(g * Lref) # Reference velocity
    Tref = Lref / Uref # Reference time

    # Geometry is fixed across the density sweep so only the reduced mass
    # changes.
    Lm = 1.0 # Membrane length
    H = Lm / 2 # Height of the domain
    domain = 1.1 * Lm # Length of the required domain
    damp = 5.0 * H # Inlet and outlet damping length
    Lf = domain + 2 * damp # Total domain length
    hs = 0.01 / Lref # Thickness of the membrane
    meshpath = WSI2DMesh1.create_mesh(ranks, H, damp, domain, Lm)

    # Damping-zone extents used by the inlet/outlet absorbing layers.
    Lfd = damp
    Lfd1 = 0.5 * H
    Ld = Lf - damp
    Ld1 = Lf - 0.5 * H

    # Time integration parameters.
    ρ∞ = 0.5
    t0 = 0.0 / Tref
    tF = 12.0 / Tref
    dt = 0.1 / Tref

    # Physical parameters. The sweep scales the reference reduced mass.
    M_standard = 0.045
    density_ratios = [0.0001, 0.001, 0.01, 0.1]
    M_sweep = density_ratios .* M_standard # Sweep over reduced mass values
    τ = 0.025 # Non dimensional pretension parameter

    # Incident-wave parameters
    η₀_dim = 0.1 # Wave amplitude [m]
    ω_dim = 2.4 # Angular frequency [rad/s]
    # Nondimensional wave parameters
    η₀ = η₀_dim / Lref
    ω = ω_dim * Tref
    # Wavenumber from dimensional finite-depth dispersion relation
    kλ_dim = find_zero(
      k -> g * k * tanh(k * H * Lref) - ω_dim^2,
      (0.01, 10.0),
    )
    kλ = kλ_dim * Lref # Nondimensional wavenumber
    ϕ = 0.0

    # Post-processing controls and in-memory sweep summary.
    vtkoutput = true
    density_sweep_data = []

    for (density_ratio, M) in zip(density_ratios, M_sweep)
      params = WSI2D_params(
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
        (reduced_mass = M, data = data),
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

    # Non dimensionalization parameters
    Lref = 20.0 # Physical Membrane length for the baseline case
    g = 9.81 # Acceleration due to gravity
    Uref = sqrt(g * Lref) # Reference velocity
    Tref = Lref / Uref # Reference time

    # Base geometry. Each membrane length gets its own mesh inside the loop.
    H = 0.5 # Height of the domain
    length_ratios = [1.0, 5.0, 10.0, 100.0]
    Lm_sweep = length_ratios # Sweep over membrane lengths
    hs = 0.01 / Lref # Thickness of the membrane

    # Time integration parameters.
    ρ∞ = 0.5
    t0 = 0.0 / Tref
    tF = 12.0 / Tref
    dt = 0.1 / Tref

    # Physical parameters held fixed across the length sweep.
    M = 0.045 # Non dimensionalized reduced mass parameter
    τ = 0.025 # Non dimensional pretension parameter

    # Incident-wave parameters. kλ is computed from the dimensional finite-depth
    # dispersion relation for the chosen frequency.
    η₀_dim = 0.1 # Wave amplitude [m]
    ω_dim = 2.4 # Angular frequency [rad/s]
    η₀ = η₀_dim / Lref
    ω = ω_dim * Tref
    kλ_dim = find_zero(
      k -> g * k * tanh(k * H * Lref) - ω_dim^2,
      (0.01, 10.0),
    )
    kλ = kλ_dim * Lref # Nondimensional wavenumber
    ϕ = 0.0

    # Post-processing controls and in-memory sweep summary.
    vtkoutput = true
    length_sweep_data = []

    for (length_ratio, Lm) in zip(length_ratios, Lm_sweep)
      damp = 5 * H # Inlet and outlet damping length

      # Keep a small clearance around the membrane while varying its length.
      domain = Lm + 0.1 # Length of the required domain
      Lf = domain + 2 * damp # Total domain length

      # Damping-zone extents for this mesh.
      Lfd = damp
      Ld = Lf - damp
      Lfd1 = 0.5 * H
      Ld1 = Lf - 0.5 * H

      meshpath = WSI2DMesh1.create_mesh(ranks, H, damp, domain, Lm)

      params = WSI2D_params(
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

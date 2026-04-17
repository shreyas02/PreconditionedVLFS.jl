module WSI3DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson, TimerOutputs
using Plots, DataFrames, Roots

include("mesh_wsi_3d_warmup.jl")
using .WSI3DMesh_warmup

include("mesh_wsi_3d_1.jl")
using .WSI3DMesh1

function warmup()
    with_mpi() do distribute
        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1)
        ranks = distribute(LinearIndices((prod(parts),)))

        # Construction of the required parameters

        # Case number
        case_name = "warmup"

        # Geometric parameters
        H = 2.0
        Lm = 1.0
        Lf = 2.0
        Ly = 2.0
        hs = 0.01
        meshpath = WSI3DMesh_warmup.create_mesh(ranks)

        # Damping parameters
        Lfd = 0.0
        Lfd1 = 0.0
        Ld = 2.0
        Ld1 = 2.0   

        # Temporal parameters
        ρ∞ = 0.5
        t0 = 0.0
        tF = 0.1
        dt = 0.1

        # Physical parameters
        ρf = 1000.0 # Fluid density
        ρs = 100.0 # Solid density
        g = 9.81 # Acceleration due to gravity
        T = 0.9 * ρf * g # Solid stiffness parameter

        # Wave parameters
        kλ = 3.0
        η₀ = 0.01
        ϕ = 0.0
        ω = sqrt(g * kλ * tanh(kλ * H)) # Wave frequency

        # Post-processing parameters
        vtkoutput = false

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

        # Loading the case and running the code
        #  not using produce_or_load here since we want to ensure the warmup runs every times
        _ = wsi3d(distribute, parts, case)
    end
end

function case_1()
    with_mpi() do distribute

        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1)
        ranks = distribute(LinearIndices((prod(parts),)))

        # Function to run the code
        function run_src(params ::WSI3D_params)
            
            solver_stats = wsi3d(distribute, parts, params)

            config_dict = Dict(string(k) => v for (k, v) in DrWatson.struct2dict(params))

            return merge(
                config_dict,
                Dict(
                    "fluid" => solver_stats.fluid,
                    "solid" => solver_stats.solid,
                    "freesurface" => solver_stats.freesurface,
                    "outer" => solver_stats.outer,
                ))
        end

        # Construction of the required parameters

        # Case number
        case_name = "case_1"

        # Geometric parameters
        H = 1.1
        Lm = 2
        Lf = 9 * pi * H
        Ly = H
        hs = 0.01

        # Damping parameters
        Lfd = 7.5 * Lm
        Lfd1 = 0.5 * Lm
        Ld = Lf - Lfd
        Ld1 = Lf - 0.5 * Lm

        # Temporal parameters
        ρ∞ = 0.5
        t0 = 0.0
        tF = 0.2
        dt = 0.1

        # Physical parameters
        ρf = 1000.0 # Fluid density
        ρs = (ρf * Lm * 0.045) / hs # Solid density
        g = 9.81 # Acceleration due to gravity
        T = 0.025 * ρf * g * Lm^2 # Solid stiffness parameter

        # Wave parameters
        η₀ = 0.1
        ω = 2 # Wave frequency
        kλ = find_zero(kλ -> g * kλ * tanh(kλ * H) - ω^2, (0.01, 10.0))
        ϕ = 0.0

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
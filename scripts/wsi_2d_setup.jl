module WSI2DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson, DataFrames
using Plots, TimerOutputs

include("mesh_wsi_2d_1.jl")
using .WSI2DMesh1

function warmup()
    with_mpi() do distribute
        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1)
        ranks = distribute(LinearIndices((prod(parts),)))

        # Construction of the required parameters

        # Case number
        case = "warmup"

        H = 1.0
        meshpath = WSI2DMesh1.create_mesh(ranks, H)

        case = WSI2D_params(
            nprocs = MPI.Comm_size(MPI.COMM_WORLD),
            rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
            case = case,
            H = H,
            meshpath = meshpath,
            ρ∞ = 0.5,
            t0 = 0.0,
            tF = 0.1,
            dt = 0.1,
            vtkoutput = false,
        )

        # Loading the case and running the code
        #  not using produce_or_load here since we want to ensure the warmup runs every times
        _ = wsi2d(distribute, parts, case)
    end
end

function case_1()
    with_mpi() do distribute

        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1)
        ranks = distribute(LinearIndices((prod(parts),)))

        # Function to run the code
        function run_src(params ::WSI2D_params)
            
            solver_stats = wsi2d(distribute, parts, params)

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
        case = "case_1"

        # Geometric parameters
        H = 10.0
        meshpath = WSI2DMesh1.create_mesh(ranks, H)
        Lf = 48*H
        
        # Temoral parameters
        ρ∞ = 0.5
        t0 = 0.0
        tF = 0.2
        dt = 0.1

        # Damping parameters 
        Lfd = 7.5*(2*H)
        Lfd1 = 0.5*(2*H)
        Ld = Lf - Lfd
        Ld1 = Lf - 0.5*(2*H)

        # Wave Parameters
        kλ = 1.0
        η₀ = 0.01

        case = WSI2D_params(
            nprocs = MPI.Comm_size(MPI.COMM_WORLD),
            rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
            case = case,
            H = H,
            meshpath = meshpath,
            Lf = Lf,
            ρ∞ = ρ∞,
            t0 = t0,
            tF = tF,
            dt = dt,
            Lfd = Lfd,
            Lfd1 = Lfd1,
            Ld = Ld,
            Ld1 = Ld1,
            kλ = kλ,
            η₀ = η₀,
            vtkoutput = true,
        )

        path = mkpath("$(datadir("wsi_2d", "case_1"))")
        filename = savename(case; ignores = [:meshpath])

        produce_or_load(run_src, case, path; filename = filename)
    end
end

end
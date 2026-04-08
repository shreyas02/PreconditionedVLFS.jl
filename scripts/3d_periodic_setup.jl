module Periodic3DSetup

using PreconditionedVLFS
using PartitionedArrays, MPI
using DrWatson, DataFrames
using Plots, TimerOutputs

function warmup()
    with_mpi() do distribute

        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1,1)
      
        # Construction of case parameters
        case = Periodic3D_params(
            iter = 0,
            case = "warmup",
            H = 1.0, # Height of the domain
            Lf = 1.0, # Length of the entire domain
            Ly = 1.0, # Width of the domain
            ρ∞ = 0.5,
            t0 = 0.0,
            tF = 0.1,
            dt = 0.1,
            vtkoutput = false,
        ) 

        # Loading the case and running the code
        #  not using produce_or_load here since we want to ensure the warmup runs every times
        _ = periodic3D(distribute, parts, case)
    end
end

function strong_scaling()
    
    with_mpi() do distribute
        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1,1)
        ranks = distribute(LinearIndices((prod(parts),)))

        # Function to run the code
        function run_src(case ::Periodic3D_params)
            solver_stats = periodic3D(distribute, parts, case)
            # Convert DrWatson's Symbol keys to String keys
            config_dict = Dict(string(k) => v for (k, v) in DrWatson.struct2dict(case))
            return merge(
                config_dict,
                Dict(
                    "fluid" => solver_stats.fluid,
                    "solid" => solver_stats.solid,
                    "outer" => solver_stats.outer,
                ))
        end
        # Strong scaling test case | Construction of the required parameters
        # We run a single case for 10 times and then take the best time to plot the strong scaling curve
        for i = 1:10
            case = Periodic3D_params(
                nprocs = MPI.Comm_size(MPI.COMM_WORLD),
                rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
                iter = i,
                case = "strong_scaling",
                H = 1.1, # Height of the domain
                Lf = 10*π, # Length of the entire domain
                Ly = 2, # Width of the domain
                dt = 0.1,
                tF = 0.2,
            )

            # Path for this test case
            path = mkpath("$(datadir("periodic3d","strong_scaling"))")
            # Run the code
            data , _ = produce_or_load(run_src,case,path)
        end
        
        ## Post processing of the data

        runtime_seconds(timer,field) = TimerOutputs.time(timer["Time Stepping Loop Aggregate"][field]) / 1e9
        path = datadir("periodic3d","strong_scaling")
        df = collect_results(path)::DataFrame

        # First Solve
        df.first_solve = [ismissing(x) ? missing : runtime_seconds(x.timer,"First solve") for x in df.outer]
        df_scaled_fs = combine(groupby(df, [:nprocs, :iter]),:first_solve => maximum => :first_solve,)
        df_scaled_fs = combine(groupby(df_scaled_fs, :nprocs),:first_solve => minimum => :first_solve,)
        sort!(df_scaled_fs, :nprocs)

        # Subsequent Solve
        df.subsequent_solve = [ismissing(x) ? missing : runtime_seconds(x.timer,"Subsequent Solves") for x in df.outer]
        df_scaled_ss = combine(groupby(df, [:nprocs, :iter]),:subsequent_solve => maximum => :subsequent_solve,)
        df_scaled_ss = combine(groupby(df_scaled_ss, :nprocs),:subsequent_solve => minimum => :subsequent_solve,)
        sort!(df_scaled_ss, :nprocs)

        # Plotting the strong scaling curve
        if i_am_main(ranks)
            # First solve
            nprocs_fs = df_scaled_fs.nprocs
            wall_fs = df_scaled_fs.first_solve
            speedup = wall_fs[1] ./ wall_fs
            ideal_speedup = nprocs_fs ./ nprocs_fs[1]

            p1 = plot(nprocs_fs, wall_fs; label = "Wall time", xlabel = "MPI Processes", ylabel = "Wall Time (s)")
            p2 = plot(nprocs_fs, speedup; label = "Measured speedup", xlabel = "MPI Processes", ylabel = "Speedup")
            plot!(p2, nprocs_fs, ideal_speedup; label = "Ideal speedup", linestyle = :dash)
            fig = plot(p1, p2; layout = (1, 2))
            savefig(fig, plotsdir("periodic3d_strong_scaling_first_solve.png"))
            
            # Subsequent solve
            nprocs_ss = df_scaled_ss.nprocs
            wall_ss = df_scaled_ss.subsequent_solve
            speedup = wall_ss[1] ./ wall_ss
            ideal_speedup = nprocs_ss ./ nprocs_ss[1]

            p3 = plot(nprocs_ss, wall_ss; label = "Wall time", xlabel = "MPI Processes", ylabel = "Wall Time (s)")
            p4 = plot(nprocs_ss, speedup; label = "Measured speedup", xlabel = "MPI Processes", ylabel = "Speedup")
            plot!(p4, nprocs_ss, ideal_speedup; label = "Ideal speedup", linestyle = :dash)
            fig = plot(p3, p4; layout = (1, 2))
            savefig(fig, plotsdir("periodic3d_strong_scaling_subsequent_solves.png"))
        end
    end
end

function weak_scaling()
    
    with_mpi() do distribute
        # Generate ranks
        parts = (MPI.Comm_size(MPI.COMM_WORLD),1,1)
        ranks = distribute(LinearIndices((prod(parts),)))

        # Function to run the code
        function run_src(case ::Periodic3D_params)
            solver_stats = periodic3D(distribute, parts, case)
            # Convert DrWatson's Symbol keys to String keys
            config_dict = Dict(string(k) => v for (k, v) in DrWatson.struct2dict(case))
            return merge(
                config_dict,
                Dict(
                    "fluid" => solver_stats.fluid,
                    "solid" => solver_stats.solid,
                    "outer" => solver_stats.outer,
                ))
        end
        # Weak scaling test case | Construction of the required parameters
        # We run a single case for 10 times and then take the best time to plot the weak scaling curve
        for i = 1:10
            case = Periodic3D_params(
                nprocs = MPI.Comm_size(MPI.COMM_WORLD),
                rank = MPI.Comm_rank(MPI.COMM_WORLD) + 1,
                case = "weak_scaling",
                iter = i,
                H = 1.1, # Height of the domain
                Lf = 10*π*(MPI.Comm_size(MPI.COMM_WORLD)+1), # Length of the entire domain
                Ly = 2, # Width of the domain
                dt = 0.1,
                tF = 0.2,
            )

            # Path for this test case
            path = mkpath("$(datadir("periodic3d","weak_scaling"))")
            # Run the code
            data , _ = produce_or_load(run_src,case,path)
        end
        
        ## Post processing of the data

        runtime_seconds(timer,field) = TimerOutputs.time(timer["Time Stepping Loop Aggregate"][field]) / 1e9
        path = datadir("periodic3d","weak_scaling")
        df = collect_results(path)::DataFrame

        # First Solve
        df.first_solve = [ismissing(x) ? missing : runtime_seconds(x.timer,"First solve") for x in df.outer]
        df_scaled_fs = combine(groupby(df, [:nprocs, :iter]),:first_solve => maximum => :first_solve,)
        df_scaled_fs = combine(groupby(df_scaled_fs, :nprocs),:first_solve => minimum => :first_solve,)
        sort!(df_scaled_fs, :nprocs)

        # Subsequent Solve
        df.subsequent_solve = [ismissing(x) ? missing : runtime_seconds(x.timer,"Subsequent Solves") for x in df.outer]
        df_scaled_ss = combine(groupby(df, [:nprocs, :iter]),:subsequent_solve => maximum => :subsequent_solve,)
        df_scaled_ss = combine(groupby(df_scaled_ss, :nprocs),:subsequent_solve => minimum => :subsequent_solve,)
        sort!(df_scaled_ss, :nprocs)

        # Plotting the weak scaling curve
        if i_am_main(ranks)
            # First solve
            nprocs_fs = df_scaled_fs.nprocs
            wall_fs = df_scaled_fs.first_solve
            ideal_wall_fs = fill(wall_fs[1], length(wall_fs))
            efficiancy_fs = wall_fs[1] ./ wall_fs

            p1 = plot(nprocs_fs, wall_fs; label = "Measured wall time", xlabel = "MPI Processes", ylabel = "Wall Time (s)")
            plot!(p1, nprocs_fs, ideal_wall_fs; label = "Ideal wall time", linestyle = :dash)
            p2 = plot(nprocs_fs, efficiancy_fs; label = "Efficiency", xlabel = "MPI Processes", ylabel = "Efficiency")
            fig = plot(p1, p2; layout = (1, 2))
            savefig(fig, plotsdir("periodic3d_weak_scaling_first_solve.png"))
            
            # Subsequent solve
            nprocs_ss = df_scaled_ss.nprocs
            wall_ss = df_scaled_ss.subsequent_solve
            ideal_wall_ss = fill(wall_ss[1], length(wall_ss))
            efficiancy_ss = wall_ss[1] ./ wall_ss

            p3 = plot(nprocs_ss, wall_ss; label = "Measured wall time", xlabel = "MPI Processes", ylabel = "Wall Time (s)")
            plot!(p3, nprocs_ss, ideal_wall_ss; label = "Ideal wall time", linestyle = :dash)
            p4 = plot(nprocs_ss, efficiancy_ss; label = "Efficiency", xlabel = "MPI Processes", ylabel = "Efficiency")
            fig = plot(p3, p4; layout = (1, 2))
            savefig(fig, plotsdir("periodic3d_weak_scaling_subsequent_solves.png"))
        end
    end
end

end # module
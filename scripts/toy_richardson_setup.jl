module ToyRichardsonSetup

using PreconditionedVLFS
using DrWatson
using Plots

function warmup()

    # Construction of the required parameters
    case = ToyRichardson_params(
        vtkoutput = false,
        tF = 0.001,
        dt = 0.001,
        length = 1.0,
        height = 1.0,
    )
    # Run the code
    #  not using produce_or_load here since we want to ensure the warmup runs every times
    _ = toyrichardson(0.5, case)
end

function comparison()

    function run_src(params ::ToyRichardson_params)
        relaxation = collect(0:0.01:1)
        iterations = Vector{Int}(undef, length(relaxation))
        for i = 1:101
            outer_num_iter, outer_iter_array, outer_residuals = toyrichardson(relaxation[i], params)
            iterations[i] = outer_num_iter - 1
        end
        return Dict(
            "relaxation_values" => relaxation,
            "iterations" => iterations)
    end

    # Path for this test case
    path = mkpath("$(datadir("toy_richardson"))")
    
    ################################
    # Comparison of robin Parameters
    _αf_opt(ρ∞, dt, ρs, hs, L) = (ρs*hs*(1-((2 * ρ∞ - 1) / (1 + ρ∞))))/(dt*(1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))*(1-(ρ∞ / (1 + ρ∞)))) + (L*dt*((1 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))^2 / 4))/(1 / 2 - (2 * ρ∞ - 1) / (1 + ρ∞) + ρ∞ / (1 + ρ∞))
    
    case1 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.001,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 45.0,
        hs = 0.1,
        L = 1e5,
        αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
        αs = 0.0) # Best case for Robin-Neumann coupling 
    
    case2 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.001,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 45.0,
        hs = 0.1,
        L = 1e5,
        αf = 1e5,
        αs = 0.0) # Dirichlet-Neumann coupling

    case3 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.001,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 45.0,
        hs = 0.1,
        L = 1e5,
        αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
        αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5)) # Robin-Robin coupling

    case4 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.001,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 45.0,
        hs = 0.1,
        L = 1e5,
        αf = 0.1 * _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
        αs = 0.0) # Sub-optimal Robin-Neumann coupling

    case5 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.001,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 45.0,
        hs = 0.1,
        L = 1e5,
        αf = 10.0 * _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
        αs = 0.0) # Sub-optimal Robin-Neumann coupling

    robincomparison_cases = [case1, case2, case3, case4, case5]

    # Setting up the plot for the comparison of the robin parameters
    p = plot(
        title = "Robin Parameter Comparison",
        xlabel = "Relaxation Parameter (ω)",
        ylabel = "Number of Iterations",
        legend = :topright, # You can change this to :topleft, :outerright, etc.
        linewidth = 2,
        ylims = (0, 49),
    )
    for case in robincomparison_cases
        data, _ = produce_or_load(run_src, case, path)
        αf = case.αf
        αs = case.αs
        iterations = data["iterations"]
        relaxation_values = data["relaxation_values"]
        legend_label = "αf = $(sci_str(αf)), αs = $(sci_str(αs))"
        plot!(p, relaxation_values, iterations, label = legend_label)
        println("Completed case with αf = $(round(αf, sigdigits=2)), αs = $(round(αs, sigdigits=2))")
    end
    save_path = plotsdir("toy_robin_comparison.png")
    savefig(p, save_path)

    ################################################################
    # Comparison of the robin-neumann method for different densities
    
    case6 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.003,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 45.0,
        hs = 0.1,
        L = 1e5,
        αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
        αs = 0.0)

    case7 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.003,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 10.0,
        hs = 0.1,
        L = 1e5,
        αf = _αf_opt(0.5, 0.001, 10.0, 0.1, 1e5),
        αs = 0.0)
    
    case8 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.003,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 1.0,
        hs = 0.1,
        L = 1e5,
        αf = _αf_opt(0.5, 0.001, 1.0, 0.1, 1e5),
        αs = 0.0)

    case9 = ToyRichardson_params(
        vtkoutput = false,
        ρ∞ = 0.5,
        tF = 0.003,
        dt = 0.001, 
        ρf = 1.0,
        ρs = 0.1,
        hs = 0.1,
        L = 1e5,
        αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
        αs = 0.0)

    densitycomparison_cases = [case6, case7, case8, case9]

    # Setting up the plot for the comparison of solver convergence for different densities
    p = plot(
        title = "Density Comparison for Robin-Neumann Coupling",
        xlabel = "Relaxation Parameter (ω)",
        ylabel = "Number of Iterations",
        legend = :topright, # You can change this to :topleft, :outerright, etc.
        linewidth = 2,
        ylims = (0, 49),
    )

    for case in densitycomparison_cases
        data, _ = produce_or_load(run_src, case, path)
        iterations = data["iterations"]
        relaxation_values = data["relaxation_values"]
        legend_label = "ρs = $(round(case.ρs, sigdigits=2)) and ρf = $(round(case.ρf, sigdigits=2))"
        plot!(p, relaxation_values, iterations, label = legend_label)
        println("Completed case with ρs = $(round(case.ρs, sigdigits=2)) and ρf = $(round(case.ρf, sigdigits=2))")
    end
    save_path = plotsdir("toy_density_comparison.png")
    savefig(p, save_path)
end
    
end # module
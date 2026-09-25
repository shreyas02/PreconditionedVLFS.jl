module ToyRichardsonSetup

using PreconditionedVLFS
using DrWatson

function comparison()

  # Case map:
  #   1-5:   robinsweep: Robin parameter sweep at fixed density and length.
  #   6-9:   robinneumann_densitysweep: Robin-Neumann density sweep.
  #   10-13: dirichletneumann_densitysweep: Dirichlet-Neumann density sweep.
  #   14-17: dirichletrobin_densitysweep: Dirichlet-Robin density sweep.
  #   18-21: robinrobin_densitysweep: Robin-Robin density sweep.
  #   23-26: robinneumann_lengthsweep: Robin-Neumann length sweep.
  #   27-30: dirichletneumann_lengthsweep: Dirichlet-Neumann length sweep.
  #   31-34: dirichletrobin_lengthsweep: Dirichlet-Robin length sweep.
  #   35-38: robinrobin_lengthsweep: Robin-Robin length sweep.
  #
  # Each case is stored through DrWatson's `produce_or_load`.  The
  # `case` field is included in the generated data filename.

  # Solver runner used by DrWatson.  Each case is evaluated for all
  # Richardson relaxation values, and only the convergence trace is stored.
  function run_src(params::ToyRichardsonParams)
    relaxation = collect(0:0.01:1)
    iterations = Vector{Int}(undef, length(relaxation))
    for i = 1:101
      outer_num_iter, _, _ = toyrichardson(relaxation[i], params)
      iterations[i] = outer_num_iter - 1
    end
    return Dict(
      "relaxation_values" => relaxation,
      "iterations" => iterations,
    )
  end

  # Small log helpers keep the experiment loops compact and consistent.
  function log_robin_case(case)
    αf = round(case.αf, sigdigits = 2)
    αs = round(case.αs, sigdigits = 2)
    println(
      "Completed case $(case.case) with " *
      "αf = $αf, αs = $αs",
    )
  end

  function log_density_case(case)
    ρs = round(case.ρs, sigdigits = 2)
    ρf = round(case.ρf, sigdigits = 2)
    println(
      "Completed case $(case.case) with " *
      "ρs = $ρs and ρf = $ρf",
    )
  end

  function log_length_case(case)
    length = round(case.length, sigdigits = 2)
    ρs = round(case.ρs, sigdigits = 2)
    ρf = round(case.ρf, sigdigits = 2)
    println(
      "Completed case $(case.case) with " *
      "length = $length, ρs = $ρs and ρf = $ρf",
    )
  end

  # Analytic optimal Robin parameter used as a baseline in the sweeps.
  function _αf_opt(ρ∞, dt, ρs, hs, L)
    γ = (2 * ρ∞ - 1) / (1 + ρ∞)
    β = ρ∞ / (1 + ρ∞)
    denom = 1 / 2 - γ + β
    return (
      (ρs * hs * (1 - γ)) / (dt * denom * (1 - β)) +
      (L * dt * ((1 - γ + β)^2 / 4)) / denom
    )
  end

  # Path for this test case
  path = mkpath("$(datadir("toy_richardson"))")

  # ------------------------------------------------------------------
  # Cases 1-5: Robin parameter comparison.
  #
  # Case 1 is the optimal Robin-Neumann baseline.  Cases 2 and 3 switch
  # coupling type, while cases 4 and 5 perturb αf around the optimum.

  case1 = ToyRichardsonParams(
    case = "toyrichardson_robinsweep_1",
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

  case2 = ToyRichardsonParams(
    case = "toyrichardson_robinsweep_2",
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

  case3 = ToyRichardsonParams(
    case = "toyrichardson_robinsweep_3",
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

  case4 = ToyRichardsonParams(
    case = "toyrichardson_robinsweep_4",
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

  case5 = ToyRichardsonParams(
    case = "toyrichardson_robinsweep_5",
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

  for case in robincomparison_cases
    produce_or_load(run_src, case, path)
    log_robin_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 6-9: Robin-Neumann density comparison.
  #
  # The coupling remains Robin-Neumann and αf is recomputed for each
  # structural density.

  case6 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_densitysweep_6",
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

  case7 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_densitysweep_7",
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

  case8 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_densitysweep_8",
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

  case9 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_densitysweep_9",
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

  for case in densitycomparison_cases
    produce_or_load(run_src, case, path)
    log_density_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 10-13: Dirichlet-Neumann density comparison.
  #
  # A large αf approximates the Dirichlet side of the coupling.

  case10 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_densitysweep_10",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  case11 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_densitysweep_11",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 10.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  case12 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_densitysweep_12",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 1.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  case13 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_densitysweep_13",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 0.1,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  densitycomparison_cases = [case10, case11, case12, case13]

  for case in densitycomparison_cases
    produce_or_load(run_src, case, path)
    log_density_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 14-17: Dirichlet-Robin density comparison.
  #
  # These cases keep the fluid side near Dirichlet and use a Robin
  # parameter on the structural side.

  case14 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_densitysweep_14",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  case15 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_densitysweep_15",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 10.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  case16 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_densitysweep_16",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 1.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  case17 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_densitysweep_17",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 0.1,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  densitycomparison_cases = [case14, case15, case16, case17]

  for case in densitycomparison_cases
    produce_or_load(run_src, case, path)
    log_density_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 18-21: Robin-Robin density comparison.
  #
  # Both sides use Robin parameters while the density ratio is varied.

  case18 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_densitysweep_18",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  case19 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_densitysweep_19",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 10.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  case20 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_densitysweep_20",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 1.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  case21 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_densitysweep_21",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    ρf = 1.0,
    ρs = 0.1,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 45.0, 0.1, 1e5))

  densitycomparison_cases = [case18, case19, case20, case21]

  for case in densitycomparison_cases
    produce_or_load(run_src, case, path)
    log_density_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 23-26: Robin-Neumann length comparison.
  #
  # The method and density are fixed while the channel length is varied.

  case23 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_lengthsweep_23",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 1.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = 0.0)

  case24 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_lengthsweep_24",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 12.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = 0.0)

  case25 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_lengthsweep_25",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 24.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = 0.0)

  case26 = ToyRichardsonParams(
    case = "toyrichardson_robinneumann_lengthsweep_26",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 48.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = 0.0)

  lengthcomparison_cases = [case23, case24, case25, case26]

  for case in lengthcomparison_cases
    produce_or_load(run_src, case, path)
    log_length_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 27-30: Dirichlet-Neumann length comparison.
  #
  # This repeats the length sweep with Dirichlet-Neumann coupling.

  case27 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_lengthsweep_27",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 6.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  case28 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_lengthsweep_28",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 12.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  case29 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_lengthsweep_29",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 24.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  case30 = ToyRichardsonParams(
    case = "toyrichardson_dirichletneumann_lengthsweep_30",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 48.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = 0.0)

  lengthcomparison_cases = [case27, case28, case29, case30]

  for case in lengthcomparison_cases
    produce_or_load(run_src, case, path)
    log_length_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 31-34: Dirichlet-Robin length comparison.
  #
  # This repeats the length sweep with Dirichlet-Robin coupling.

  case31 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_lengthsweep_31",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 6.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  case32 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_lengthsweep_32",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 12.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  case33 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_lengthsweep_33",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 24.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  case34 = ToyRichardsonParams(
    case = "toyrichardson_dirichletrobin_lengthsweep_34",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 48.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = 1e7,
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  lengthcomparison_cases = [case31, case32, case33, case34]

  for case in lengthcomparison_cases
    produce_or_load(run_src, case, path)
    log_length_case(case)
  end

  # ------------------------------------------------------------------
  # Cases 35-38: Robin-Robin length comparison.
  #
  # This repeats the length sweep with Robin parameters on both sides.

  case35 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_lengthsweep_35",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 6.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  case36 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_lengthsweep_36",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 12.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  case37 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_lengthsweep_37",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 24.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  case38 = ToyRichardsonParams(
    case = "toyrichardson_robinrobin_lengthsweep_38",
    vtkoutput = false,
    ρ∞ = 0.5,
    tF = 0.003,
    dt = 0.001,
    length = 48.0,
    ρf = 1.0,
    ρs = 45.0,
    hs = 0.1,
    L = 1e5,
    αf = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5),
    αs = _αf_opt(0.5, 0.001, 0.1, 0.1, 1e5))

  lengthcomparison_cases = [case35, case36, case37, case38]

  for case in lengthcomparison_cases
    produce_or_load(run_src, case, path)
    log_length_case(case)
  end
end

end # module

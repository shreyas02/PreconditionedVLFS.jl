using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("toy_richardson_setup.jl"))
using .ToyRichardsonSetup

isempty(ARGS) && error("Pass one of: comparison, all")

case = ARGS[1]
if case == "comparison"
  ToyRichardsonSetup.comparison()
elseif case == "all"
  ToyRichardsonSetup.comparison()
else
  println(stderr, "Unknown Toy Richardson case: $case")
  exit(1)
end

using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("2d_periodic_setup.jl"))
using .Periodic2DSetup

isempty(ARGS) && error("Pass one of: strong_scaling, weak_scaling, all")

case = ARGS[1]
if case == "strong_scaling"
  Periodic2DSetup.strong_scaling()
elseif case == "weak_scaling"
  Periodic2DSetup.weak_scaling()
elseif case == "all"
  Periodic2DSetup.strong_scaling()
  Periodic2DSetup.weak_scaling()
else
  println(stderr, "Unknown Periodic 2D case: $case")
  exit(1)
end

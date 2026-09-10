using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("3d_periodic_setup.jl"))
using .Periodic3DSetup

isempty(ARGS) && error("Pass one of: strong_scaling, weak_scaling, all")

case = ARGS[1]
if case == "strong_scaling"
  Periodic3DSetup.strong_scaling()
elseif case == "weak_scaling"
  Periodic3DSetup.weak_scaling()
elseif case == "all"
  Periodic3DSetup.strong_scaling()
  Periodic3DSetup.weak_scaling()
else
  println(stderr, "Unknown Periodic 3D case: $case")
  exit(1)
end

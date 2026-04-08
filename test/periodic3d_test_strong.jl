using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("3d_periodic_setup.jl"))
using .Periodic3DSetup

Periodic3DSetup.strong_scaling()
using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("2d_periodic_setup.jl"))
using .Periodic2DSetup

Periodic2DSetup.strong_scaling()
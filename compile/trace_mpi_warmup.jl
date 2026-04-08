using DrWatson
@quickactivate "PreconditionedVLFS"

using PreconditionedVLFS

include(scriptsdir("wsi_2d_setup.jl"))
include(scriptsdir("wsi_3d_setup.jl"))
include(scriptsdir("2d_periodic_setup.jl"))
include(scriptsdir("3d_periodic_setup.jl"))

using .WSI2DSetup
using .WSI3DSetup
using .Periodic2DSetup
using .Periodic3DSetup

WSI2DSetup.warmup()
WSI3DSetup.warmup()
Periodic2DSetup.warmup()
Periodic3DSetup.warmup()

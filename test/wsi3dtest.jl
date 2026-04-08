using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("wsi_3d_setup.jl"))
using .WSI3DSetup

WSI3DSetup.case_1()
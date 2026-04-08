using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("wsi_2d_setup.jl"))
using .WSI2DSetup

WSI2DSetup.case_1()
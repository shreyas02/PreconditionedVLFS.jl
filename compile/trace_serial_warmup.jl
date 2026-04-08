using DrWatson
@quickactivate "PreconditionedVLFS"

using PreconditionedVLFS

include(scriptsdir("toy_richardson_setup.jl"))
using .ToyRichardsonSetup

ToyRichardsonSetup.warmup()

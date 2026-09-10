using DrWatson
using PreconditionedVLFS

include(scriptsdir("wsi_2d_setup.jl"))
using .WSI2DSetup

isempty(ARGS) && error("Pass one of: case_1, length_sweep, density_sweep, all")

case = ARGS[1]
if case == "case_1"
  WSI2DSetup.case_1()
elseif case == "length_sweep"
  WSI2DSetup.case_1_length_sweep()
elseif case == "density_sweep"
  WSI2DSetup.case_1_density_sweep()
elseif case == "all"
  WSI2DSetup.case_1()
  WSI2DSetup.case_1_length_sweep()
  WSI2DSetup.case_1_density_sweep()
else
  println(stderr, "Unknown WSI 2D case: $case")
  exit(1)
end

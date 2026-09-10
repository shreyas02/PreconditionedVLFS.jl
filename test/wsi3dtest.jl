using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("wsi_3d_setup.jl"))
using .WSI3DSetup

isempty(ARGS) && error("Pass one of: case_1, all")

case = ARGS[1]
if case == "case_1"
  WSI3DSetup.case_1()
elseif case == "all"
  WSI3DSetup.case_1()
else
  println(stderr, "Unknown WSI 3D case: $case")
  exit(1)
end

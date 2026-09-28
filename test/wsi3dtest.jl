using DrWatson
@quickactivate "PreconditionedVLFS"
using PreconditionedVLFS

include(scriptsdir("wsi_3d_setup.jl"))
using .WSI3DSetup

isempty(ARGS) && error("Pass one of: case_1, case_2, all")

case = ARGS[1]
if case == "case_1"
  WSI3DSetup.case_1()
elseif case == "case_2"
  WSI3DSetup.case_2()
elseif case == "all"
  WSI3DSetup.case_1()
  WSI3DSetup.case_2()
else
  println(stderr, "Unknown WSI 3D case: $case")
  exit(1)
end

using MPIPreferences
using PackageCompiler

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const PRECOMPILE_STATEMENTS = joinpath(@__DIR__, "PreconditionedVLFS_precompile_statements.jl")
const SYSIMAGE_PATH = joinpath(@__DIR__, "PreconditionedVLFS.so")

if !isfile(PRECOMPILE_STATEMENTS)
    error(
        "Missing precompile statements file: $(PRECOMPILE_STATEMENTS).\n" *
        "Run `bash run/build_sysimage.sh` first to generate it.",
    )
end

const MPI_BINARY = string(MPIPreferences.binary)
const MPI_ABI = string(MPIPreferences.abi)

println("Building sysimage with MPI binary=$(MPI_BINARY), abi=$(MPI_ABI)")

mkpath(dirname(SYSIMAGE_PATH))

create_sysimage(
    [:PreconditionedVLFS];
    project = PROJECT_ROOT,
    sysimage_path = SYSIMAGE_PATH,
    precompile_statements_file = PRECOMPILE_STATEMENTS,
    incremental = false,
    include_transitive_dependencies=false,
)

println("Created sysimage at $(SYSIMAGE_PATH)")
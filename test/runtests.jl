using NeighborJoining, Test, SafeTestsets

const GROUP = get(ENV, "GROUP", "All")
const is_CI = haskey(ENV, "CI")
const KNOWN_GROUPS = ["All", "Aqua", "RNJ", "FNJ"]

GROUP in KNOWN_GROUPS || AssertionError("Env. Var. 'GROUP': $GROUP, not in testing groups\nTesting groups are $KNOWN_GROUPS")
@info "GROUP is " GROUP
@time begin

    if GROUP == "All" || GROUP == "RNJ"
        @safetestset "RegularNeighborJoining" begin
            include("testRegularNeighborJoining.jl")
        end
    end

    if GROUP == "All" || GROUP == "FNJ"
        @safetestset "FastNeighborJoining" begin
            include("testFastNeighborJoining.jl")
        end
    end

    if GROUP == "All" || GROUP == "Aqua"
        @safetestset "Aqua Q/A" begin
            using Aqua, NeighborJoining
            Aqua.find_persistent_tasks_deps(NeighborJoining)
            Aqua.test_ambiguities(NeighborJoining, recursive = false)
            Aqua.test_deps_compat(NeighborJoining)
            Aqua.test_project_extras(NeighborJoining)
            Aqua.test_stale_deps(NeighborJoining)
            Aqua.test_unbound_args(NeighborJoining)
            Aqua.test_undefined_exports(NeighborJoining)
        end
    end
end

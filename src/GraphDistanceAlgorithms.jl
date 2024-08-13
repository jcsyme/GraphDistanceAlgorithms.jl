"""
Graph distance algorithms for integration into other algorithms. The 
    GraphDistanceAlgorithms package supports the use of passing vectors to
    distance algorithms to reduce memory pressure.
    
    * The Dijkstra implementation includes both KaryHeaps (from IterativeHeaps)
        and QuickHeaps. KaryHeaps support the use of passing vectors.
"""
module GraphDistanceAlgorithms

using Distributed
using DistributedArrays
using Graphs
using SharedArrays

# custom modules TEMPORARY
using Pkg
Pkg.develop(path = "/Users/jsyme/Documents/Projects/git_jbus/IterativeHeaps.jl")
using IterativeHeaps


export bellman_ford,
       bellman_ford!,
       check_arrays_for_algorithm,
       dijkstra!,
       dijkstra_kary,
       dijkstra_kary!,
       dijkstra_quickheaps,
       select_algorithm,
       spawn_arrays



# some base types
VMOrNoth{T} = Union{Vector{T}, Matrix{T}, Nothing}


# include
include("BellmanFord.jl")
include("Dijkstra.jl")
include("Utilities.jl")


end 

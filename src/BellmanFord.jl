"""
Find graph distances from source vertices to all others using the Bellman-Ford 
    shortest paths using Kary Heaps (generic formulation of Binary heaps to 
    _k_), which allows for passing of vectors to reduce memory overhead.

# Constructs

```
bellman_ford(
    graph::AbstractGraph,
    sources::Union{Vector{Int64}, Int64};
    active::VMOrNoth{Bool} = nothing,
    dists::VMOrNoth{Int64} = nothing,
    distmx::AbstractMatrix{T} = Graphs.weights(graph),
    new_active::VMOrNoth{Bool} = nothing,
    parents::VMOrNoth{Int64} = nothing,
)
```

```
bellman_ford!(
    dists::Union{Vector{T}, Matrix{T}},
    graph::AbstractGraph,
    sources::Union{Vector{Int64}, Int64};
    active::VMOrNoth{Bool} = nothing,
    dists::VMOrNoth{Int64} = nothing,
    distmx::AbstractMatrix{T} = Graphs.weights(graph),
    new_active::VMOrNoth{Bool} = nothing,
    parents::VMOrNoth{Int64} = nothing,
)
```

##  Function Arguments

- `dists`: vector of distances to modify (keyword argument in non-inline)
- `g`: Graph to run algorithm on
- `srcs`: source vertices to use to calculate distances from
- `dict_shared_arrays`: for implementation of dijkstra on workers using
    SharedArray to perform calculations. Must include the following keys 
    (UNSAFE, no checking):
    * dists
    * heap_data
    * heap_index
    * heap_index_lookup
    * parents
    * size

##  Keyword Arguments

- `active`: Optional vector to pass for "active" in Bellman-Ford
- `distmx`: Optional distance matrix to pass. Defaults to Graph edge weights
- `heap_data`: Optional vector to use to store heap data. Must be of type T
    (type of edge weights)
- `heap_index`: Optional vector to use to store heap index
- `heap_index_lookup`: Optional vector to use to store heap lookup index
- `k`: size of K-ary heap. Defaults to 4
- `parents`: optional vector of parents to pass.

"""
function bellman_ford!(
    dists::Union{Vector{T}, Matrix{T}},
    graph::AbstractGraph,
    sources::Union{Vector{Int64}, Int64};
    active::VMOrNoth{Bool} = nothing,
    distmx::AbstractMatrix{T} = Graphs.weights(graph),
    new_active::VMOrNoth{Bool} = nothing,
    parents::VMOrNoth{Int64} = nothing,
) where {T<:Real}
    
    isa(sources, Int64) && (sources = [sources]);
    nvg = nv(graph)


    ##  DERIVED FROM GRAPHS.JL IMPLEMENTATION
    
    # active vertex
    active = IterativeHeaps.initialize_vector(
        active, 
        nvg, 
        Bool; 
        fill_func = zero, 
    )
    
    # new active vertex
    new_active = IterativeHeaps.initialize_vector(
        new_active, 
        nvg, 
        Bool; 
        fill_func = zero, 
    )
    
    # distances
    dists = IterativeHeaps.initialize_vector(
        dists, 
        nvg, 
        T; 
        fill_func = typemax, 
    )

    # parents
    parents = IterativeHeaps.initialize_vector(
        parents, 
        nvg, 
        Int64;
    )

    # update some source values
    active[sources] .= true
    dists[sources] .= 0
    no_changes = false

    for i in vertices(graph)
        no_changes = true
        new_active .= false

        for u in vertices(graph)
            active[u] || continue
            for v in outneighbors(graph, u)
                relax_dist = distmx[u, v] + dists[u]
                if dists[v] > relax_dist
                    dists[v] = relax_dist
                    parents[v] = u
                    no_changes = false
                    new_active[v] = true
                end
            end
        end

        no_changes && break
        active, new_active = new_active, active
    end
    
    no_changes || throw(Graphs.NegativeCycleError())
    
end

function bellman_ford(
    graph::AbstractGraph,
    sources::Union{Vector{Int64}, Int64};
    active::VMOrNoth{Bool} = nothing,
    dists::VMOrNoth{Bool} = nothing,
    distmx::AbstractMatrix{T} = Graphs.weights(graph),
    new_active::VMOrNoth{Bool} = nothing,
    parents::VMOrNoth{Int64} = nothing,
) where {T<:Real}

    nvg = nv(graph)
    
    ##  DERIVED FROM GRAPHS.JL IMPLEMENTATION
    
    # initialize distances
    dists = IterativeHeaps.initialize_vector(
        dists, 
        nvg, 
        T; 
        fill_func = typemax, 
    )

    # call bellman
    bellman_ford!(
        dists,
        graph,
        sources;
        active = active,
        distmx = distmx,
        new_active = new_active,
        parents = parents,
    )

    return dists
    
end




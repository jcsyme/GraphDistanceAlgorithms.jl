
@inline function dijkstra!(
    dists::Union{Vector{T}, Matrix{T}},
    g::AbstractGraph,
    srcs::Union{U, Vector{U}};
    heap_type::Symbol = :quickheaps,
    kwargs...
) where {T <: Real} where {U<:Integer}
    
    if heap_type == :quickheaps
        # need to set this up
        @info("dijkstra! for quickheaps currently undefined")
        #out = dijkstra_quickheaps(g, srcs; kwargs...)

    elseif heap_type == :kary
        out = dijkstra_kary!(
            dists,
            g, 
            srcs; 
            kwargs...
        )
    end

end



"""
Implement a version of Dijkstra shortest paths using Kary Heaps (generic
    formulation of Binary heaps to _k_), which allows for passing of vectors
    to reduce memory overhead.

# Constructs

```
dijkstra_kary!(
    dists::Union{Vector{T}, Matrix{T}},
    g::AbstractGraph,
    srcs::Union{U, Vector{U}};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    heap_data::VMOrNoth{T} = nothing,
    heap_index::VMOrNoth{Int64} = nothing,
    heap_index_lookup::VMOrNoth{Int64} = nothing,
    k::Int64 = 4,
    parents::VMOrNoth{Int64} = nothing,
)
```

```
function dijkstra_kary!(
    g::AbstractGraph,
    srcs::Union{U, Vector{U}},
    dict_shared_arrays::Dict{Symbol, SharedArray};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    k::Int64 = 4,
)
```

```
dijkstra_kary(
    g::AbstractGraph,
    srcs::Union{U, Vector{U}};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    dists::VMOrNoth{T} = nothing,
    heap_data::VMOrNoth{T} = nothing,
    heap_index::VMOrNoth{Int64} = nothing,
    heap_index_lookup::VMOrNoth{Int64} = nothing,
    k::Int64 = 4,
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

- `distmx`: Optional distance matrix to pass. Defaults to Graph edge weights
- `heap_data`: Optional vector to use to store heap data. Must be of type T
    (type of edge weights)
- `heap_index`: Optional vector to use to store heap index
- `heap_index_lookup`: Optional vector to use to store heap lookup index
- `k`: size of K-ary heap. Defaults to 4
- `parents`: optional vector of parents to pass.

"""
@inline function dijkstra_kary!(
    dists::Union{Vector{T}, Matrix{T}},
    g::AbstractGraph,
    srcs::Union{U, Vector{U}};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    heap_data::VMOrNoth{T} = nothing,
    heap_index::VMOrNoth{Int64} = nothing,
    heap_index_lookup::VMOrNoth{Int64} = nothing,
    k::Int64 = 4,
    parents::VMOrNoth{Int64} = nothing,
) where {T<:Real} where {U<:Integer}
    
    isa(srcs, Int64) && (srcs = [srcs]);

    # initialize vectors
    nvg = nv(g)

    dists = IterativeHeaps.initialize_vector(
        dists, 
        nvg, 
        T; 
        fill_func = typemax, 
    )
    
    parents = IterativeHeaps.initialize_vector(
        parents, 
        nvg, 
        U;
    )

    # get the heap
    H = IterativeHeaps.KaryHeap{T}(
        nvg, 
        k;
        data = heap_data,
        index = heap_index,
        index_lookup = heap_index_lookup,
    )

    for src in srcs
        @inbounds dists[src] = zero(T)
        IterativeHeaps.heap_push!(H, zero(T), src)
    end

    
    @inbounds begin
        while !Base.isempty(H)#H.size[1] > 0

            u, value = IterativeHeaps.heap_pop!(H)

            (value > dists[u]) && continue
            
            for v in outneighbors(g, u)
                
                alt = dists[u] + distmx[u, v]
                (alt >= dists[v]) && continue
                dists[v] = alt
                
                parents[v] = u
                IterativeHeaps.heap_push!(H, alt, v)

            end
            
        end
    end

    return dists
end


@inline function dijkstra_kary!(
    g::AbstractGraph,
    srcs::Union{U, Vector{U}},
    dict_shared_arrays::Dict{Symbol, SharedArray};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    k::Int64 = 4,
) where {T<:Real} where {U<:Integer}
    
    # some checks
    (myid() == 1) && (return nothing)
    isa(srcs, Int64) && (srcs = [srcs]);

    # initialize vectors
    nvg = nv(g)

    dists = dict_shared_arrays[:dists]
    kd = dists.pidx
    fill_shared_array!(
        dists;
        fill_func = typemax, 
    )
    
    parents = dict_shared_arrays[:parents]
    kp = parents.pidx
    fill_shared_array!(
        parents;
        fill_func = zero, 
    )

    sizes = dict_shared_arrays[:size]
    ks = sizes.pidx
    fill_shared_array!(
        sizes;
        fill_func = zero, 
    )

    # get the heap
    H = IterativeHeaps.KaryHeapShared{T}(
        nvg, 
        k,
        dict_shared_arrays[:heap_data],
        dict_shared_arrays[:heap_index],
        dict_shared_arrays[:heap_index_lookup],
        dict_shared_arrays[:size],
    )

    
    for src in srcs
        @inbounds dists.s[src, kd] = zero(T)
        IterativeHeaps.heap_push!(H, zero(T), src)
    end

    
    @inbounds begin
        while !Base.isempty(H)#H.size[1] > 0

            u, value = IterativeHeaps.heap_pop!(H)

            (value > dists.s[u, kd]) && continue
            
            for v in outneighbors(g, u)
                
                alt = dists.s[u, kd] + distmx[u, v]
                (alt >= dists[v, kd]) && continue
                dists.s[v, kd] = alt
                
                parents.s[v, kp] = u
                IterativeHeaps.heap_push!(H, alt, v)

            end
            
        end
    end
end



@inline function dijkstra_kary(
    g::AbstractGraph,
    srcs::Union{U, Vector{U}};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    dists::VMOrNoth{T} = nothing,
    heap_data::VMOrNoth{T} = nothing,
    heap_index::VMOrNoth{Int64} = nothing,
    heap_index_lookup::VMOrNoth{Int64} = nothing,
    k::Int64 = 4,
    parents::VMOrNoth{Int64} = nothing,
) where {T<:Real} where {U<:Integer}
    
    isa(srcs, Int64) && (srcs = [srcs]);

    # initialize vectors
    nvg = nv(g)

    dists = IterativeHeaps.initialize_vector(
        dists, 
        nvg, 
        T; 
        fill_func = typemax, 
    )
    
    # 
    dijkstra_kary!(
        dists,
        g,
        srcs;
        distmx = distmx,
        heap_data = heap_data,
        heap_index = heap_index,
        heap_index_lookup = heap_index_lookup,
        k = k,
        parents = parents,
    ) 

    return dists
end



@inline function dijkstra_quickheaps(
    g::AbstractGraph,
    srcs::Union{U, Vector{U}};
    distmx::AbstractMatrix{T} = Graphs.weights(g),
    dists::Union{Vector{Int64}, Matrix{Int64}, Nothing} = nothing,
    parents::Union{Vector{Int64}, Matrix{Int64}, Nothing} = nothing,
) where {T<:Real} where {U<:Integer}
    
    
    isa(srcs, Int64) && (srcs = [srcs]);
    
    nvg = nv(g)
    dists = IterativeHeaps.initialize_vector(
        dists, 
        nvg, 
        Int64; 
        fill_func = typemax, 
    )
    parents = IterativeHeaps.initialize_vector(
        parents, 
        nvg, 
        U; 
    )

    H = FastPriorityQueue{T}(nvg)

    for src in srcs
        @inbounds dists[src] = 0
        H[src] = zero(T)
    end
    

    while !isempty(H)

        u, value = dequeue_pair!(H)
        # u, value = FibonacciHeaps.heap_pop!(H)
        (value > dists[u]) && continue
        
        for v in outneighbors(g, u)
            @inbounds alt = dists[u] + distmx[u, v]
            (@inbounds alt >= dists[v]) && continue
            
            @inbounds dists[v] = alt
            H[v] = alt
        end
        
    end


    return dists#DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end



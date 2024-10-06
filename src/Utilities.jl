"""
Check whether a dictionary defined for arrays is specified properly.

# Constructs

```
check_arrays_for_algorithm(
    algorithm::Symbol,
    dict_arrays::Union{Dict, Nothing},
)
```


##  Function Arguments

- `algorithm`: algorithm to test
- `dict_arrays`: dictionary to pass containing cache arrays for the distance 
    algorithms
"""
function check_arrays_for_algorithm(
    algorithm::Symbol,
    dict_arrays::Union{Dict, Nothing},
)

    # check and initialize
    !(algorithm in [:bellman_ford, :dijkstra_kary, :dijkstra_quickheaps]) && (return nothing);

    if isa(dict_arrays, Nothing)
        return false
    end

    # 
    if algorithm == :bellman_ford
        keys_req = [
            :active,
            :dists,
            :new_active,
            :parents 
        ]

    elseif algorithm == :dijkstra_kary
        keys_req = [
            :dists,
            :heap_data,
            :heap_index,
            :heap_index_lookup,
            :parents,
            :size
        ]

    elseif algorithm == :dijkstra_quickheaps
        keys_req = [
            :dists,
            :parents 
        ]

    end

    out = issubset(Set(keys_req), Set(keys(dict_arrays)))

    return out
end



"""
Select an algorithm based on dimension or input. Returns a symbol

```
select_algorithm(
    algorithm::Symbol,
    n_vertices::Int64,
)
```

##  Function Arguments

- `algorithm`: algorithm to use


##  Keyword Arguments

- `iterating`: are you calling this to support an iteration? Controls the
    Dijkstra heap call.
- `n_vertices`: optional number of vertices to pass. If n_vertices < 50, and
    algorithm == :auto, selects :dijkstra
"""
function select_algorithm(
    algorithm::Symbol;
    iterating::Bool = true,
    n_vertices::Union{Int64, Nothing} = nothing,
)::Symbol

    # set algorithm
    algorithm = !(
        algorithm in [
            :auto, 
            :bellman_ford, 
            :dijkstra_kary,
            :dijkstra_quickheaps,
        ]
    ) ? :auto : algorithm

    # set the Dijkstra heap call
    algorithm_dijkstra = iterating ? :dijkstra_kary : :dijkstra_quickheaps

    if algorithm == :auto
        algorithm = (
            !isa(n_vertices, Nothing)
            ? ((n_vertices < 50) ? :dijkstra : :bellman_ford)
            : :bellman_ford
        )
    end

    alg_func = (
        algorithm == :bellman_ford
        ? :bellman_ford
        : algorithm_dijkstra
    )

    return alg_func
end



"""
Spawn Arrays, DistributedArrays, or SharedArrays for use in iterating 
    algorithms. If nworkers() == 1, spawns vectors. Otherwise, each array has
    dims (n_vertices, n_workers)

# Constructs

```
spawn_arrays(
    graph::AbstractGraph,
    dict_key_to_type::Dict{Symbol, Type};
    dict_key_to_nrow::Union{Dict{Symbol, Int64}, Nothing} = nothing,
    type::Symbol = :DistributedArray,
)
```

```
spawn_arrays(
    graph::AbstractGraph,
    algorithm::Symbol;
    type::Symbol = :DistributedArray,
)
```

##  Function Arguments

- `graph`: Graph to use for dimensions
- `algorithm`: one of the following algorithm types:
    * :bellman_ford
    * :dijkstra_kary
    * :dijkstra_quickheaps
- `dict_key_to_type`: dictionary specifying keys to type to specify. Keys 
    represent a different vector (parallelized, creating a nv x nw array) and
    types represent the data type for that array


##  Keyword Arguments

- `dict_key_to_nrow`: optional dictionary mapping a key in `dict_key_to_type` to
    a number of rows. If the key is not specified, each array defaults to 
    nv(graph) rows
- `type`: Either `:DistributeArray`, `:SharedArray`, or `:Vector`

"""
function spawn_arrays(
    graph::AbstractGraph,
    dict_key_to_type::Dict{Symbol, Type};
    dict_key_to_nrow::Union{Dict{Symbol, Int64}, Nothing} = nothing,
    type::Symbol = :DistributedArray,
)
    nvg = nv(graph)
    nw = nworkers()

    # determine output type?
    if (nw == 1) | (type == :Vector)
        dict_out = Dict{Symbol, Vector}()
        type = :Vector
        vec_out_q = true

    else
        !(type in [:DistributedArray, :SharedArray]) && (type = :DistributedArray)
        T = (type == :DistributedArray) ? DArray : SharedArray
        dict_out = Dict{Symbol, T}()
        vec_out_q = false
    end

    # iterate to build outputs
    for pair in dict_key_to_type

        # get key, type, and number of rows
        key = pair[1]
        tp = pair[2]
        m = isa(dict_key_to_nrow, Dict) ? get(dict_key_to_nrow, key, nvg) : nvg

        # build the arrays
        array = vec_out_q ? zeros(tp, m) : nothing

        # handle distributed and sharedarrays case
        if isa(array, Nothing)
            
            array = (
                (type == :DistributedArray)
                ? dzeros(tp, (m, nw), workers(), [1, nw])
                : SharedArray{tp}((m, nw); pids = workers())
            )
        end
        
        dict_out[key] = array
    end

    return dict_out
end

function spawn_arrays(
    graph::AbstractGraph,
    algorithm::Symbol;
    type::Symbol = :DistributedArray,
)
    # check and initialize
    !(algorithm in [:bellman_ford, :dijkstra_kary, :dijkstra_quickheaps]) && (return nothing);
    dict_spec = Dict{Symbol, Type}()
    dict_nrows = Dict{Symbol, Int64}()

    # all include distances and parents
    dict_spec[:dists] = Int64
    dict_spec[:parents] = Int64

    # algorithm-specific additions
    if algorithm == :bellman_ford
        dict_spec[:active] = Bool
        dict_spec[:new_active] = Bool

    elseif algorithm == :dijkstra_kary
        # add heap vectors
        dict_spec[:heap_data] = Int64
        dict_spec[:heap_index] = Int64
        dict_spec[:heap_index_lookup] = Int64
        dict_spec[:size] = Int64

        # size only needs one row
        dict_nrows[:size] = 1

    elseif algorithm == :dijkstra_quickheaps
        nothing

    end
    

    # call generic implementation
    dict_out = spawn_arrays(
        graph,
        dict_spec;
        dict_key_to_nrow = dict_nrows,
        type = type,
    )

    return dict_out
end


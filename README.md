# GraphDistanceAlgorithms.jl

## Introduction

`GraphDistanceAlgorithms.jl` includes implementations of graph distance algorithms designed to facilitate iterative parallelization using DistributedArrays to reduce memory pressure and increase speed. Currently, the package includes the followng algorithms:

- Dijkstra
- Bellman-Ford

The implementations of these algorithms are derived from code included in the `Graphs.jl` (https://juliagraphs.org/Graphs.jl/v1.5/).

The algorithms make use of `DistributedArrays` and `IterativeHeaps.jl` to pass cache arrays to heaps and other storage vectors when running in parellel; this use of preallocation is a standard Julia practice that eliminates the need to assign new memory addresses and run garbage collection during iteration. The algorithms use asynchrnous data parallelization instead of multithreading by passing each source vertex to a unique computational core.


## Use

The algorithms operate on a AbstractGraph object (can be directed or undirected graph). In general, algorithms can be run either with or without preallocated arrays; if run without, they operate in the same way that standard implementations run (e.g., see ?dijkstra, ?dijkstra!, ?bellman_ford, or ?bellman_ford! for information on input arguments).

To instantiante the algorithm for iteration, first a dictionary of input ararys must be called. These input arrays support `DistributedArrays.Darray` or `SharedArray.SharedArray` parallel array objects, though `DistributedArrays.Darray` are strongly recommended for speed. The user, however, does not need to manually instantiate these arrays. Instead, the use can build the required set of arrays in a dictionary as a function of the graph, the algorithm that will be run, and the type of parallel array to use. For example,

```
dict_arrays = spawn_arrays(
    graph::AbstractGraph,
    algorithm_run::Symbol;
    type::Symbol = :DistributedArray,
)
```

where `algorithm_run` is one of the following. 

- :bellman_ford
- :dijkstra_kary (using IterativeHeaps.jl K-Ary heap) 
- :dijkstra_quickheaps (using `QuickHeaps.jl` implentation: **NOTE: Under development**)


Note that arrays can be spawned for a single process, in which case the dictionary maps to vectors and arrays. If this occurs, then the process index [:L] (see below) does not need to be included; this only is valid for `DArray` types. See [DistributedArrays](https://juliaparallel.org/DistributedArrays.jl/stable/) for more.


###  Bellman-Ford Algorithm

To run the Bellman-Ford algorithm on graph `graph`, use two steps. First, build the arrays:

```
dict_arrays = spawn_arrays(
    graph,
    :bellman_ford;
    :DistributedArray,
)
```

Then, call the algorithm using the following inputs (for source vertex `i`) if using multiple processes:

```
bellman_ford!(
    dict_arrays[:dists][:L],
    graph, 
    i; 
    active = dict_arrays[:active][:L],
    new_active = dict_arrays[:new_active][:L],
    parents = dict_arrays[:parents][:L],
)
```


###  Dijkstra's Algorithm

To run Dijkstra's algorithm on graph `graph` using two steps. First, build the arrays for the K-Ary implementation:

```
dict_arrays = spawn_arrays(
    graph,
    :dijkstra_kary;
    :DistributedArray,
)
```

Then, call the algorithm using the following inputs (for source vertex `i`) on a distributed process (note that [:L] is the index of the node from which it is called):

```
dijkstra_kary!(
    dict_arrays[:dists][:L],
    graph, 
    i; 
    parents = dict_arrays[:parents][:L],
    heap_data = dict_arrays[:heap_data][:L],
    heap_index = dict_arrays[:heap_index][:L],
    heap_index_lookup = dict_arrays[:heap_index_lookup][:L],

)
```


See `graph_distance_algorithms_example.ipynb` for an example and comparison of benchmarks with a serial implementation.


## Data

The code only needs a Graph to work off of. This can be loaded in a Julia session using `Graphs.jl`. If using `DiscreteGraphAlgorithms.jl`, then the `graph` property of a `GraphWrapper` can be used.


## Project information

ADD HERE


## References/Bibliography

Fairbanks, James and Besan{\c{c}}on, Mathieu and Simon, Sch{\"o}lly and Hoffiman, J{\'u}lio and Eubank, Nick and Karpinski, Stefan. 2021. _JuliaGraphs/Graphs.jl: an optimized graphs package for the Julia programming language_. https://github.com/JuliaGraphs/Graphs.jl/

 

## Copyright and License

Copyright (C) <2024> RAND Corporation. This code is made available under the MIT license.

 

## Authors and Reference

James Syme, 2024.

@misc{GDA2024,
  author       = {Syme, James},
  title        = {GraphDistanceAlgorithms.jl: Graph shortest path algorithms for iteration.},
  year         = 2024,
  url = {URLHERE}
}

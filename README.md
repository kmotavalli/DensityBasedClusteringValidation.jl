# DBCV

<!-- Title -->
<h1 align="center">
DensityBasedClusteringValidation.jl
</h1>

<!-- description -->
<p align="center">
  <strong>Density-Based Clustering Validation index implementation for Julia</strong>
</p>


[![Build Status](https://github.com/kmotavalli/DensityBasedClusteringValidation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kmotavalli/DensityBasedClusteringValidation.jl/actions/workflows/CI.yml?query=branch%3Amain)

DensityBasedClusteringValidation.jl is a Julia package implementing the DBCV metric, used to map the clustering of a dataset (as obtained, for example, via the DBSCAN clustering algorithm) to an index between -1, for poorly formed clusters or wrong a points-to-cluster assignment, to +1, if presented with an obtimal clusterization of the dataset.

It identifies clusters via density variations in relation to the between-clusters density.

DBCV was proposed in:
>*Density-Based Clustering Validation. Davoud Moulavi, Pablo A. Jaskowiak, Ricardo J. G. B. Campello, Arthur Zimek, and Jörg Sander. Proceedings of the 2014 SIAM International Conference on Data Mining (SDM). 2014, 839-847*

downloadable at: [doi.org/10.1137/1.9781611973440.96]( https://doi.org/10.1137/1.9781611973440.96)

It yelds adequate results, compared to other indexes, also when evaluating data with nested clusters, concave shaped clusters, spirals, and data that is not clusterized in quadrants, while not being limited to these shapes. This is better explained in 

>*Chicco D, Sabino G, Oneto L, Jurman G. 2025. The DBCV index is more
informative than DCSI, CDbw, and VIASCKDE indices for unsupervised
clustering internal assessment of concave-shaped and density-based
clusters. PeerJ Computer Science 11:e3095*

downloadable at: [doi.org/10.7717/peerj-cs.3095](https://doi.org/10.7717/peerj-cs.3095)

## Compatibility in results with other implementations

This package for Julia tries to be similar in usage to [FelSiq/DBCV for Python](https://github.com/FelSiq/DBCV) supporting similar options, with difference better documented below, and more importantly, tries to match the calculated index value with the one derived by FelSiq/DBCV, when evaluating the same dataset and classification.

In tests, the maximum divergence between Dbcv.jl and FelSiq/DBCV on the same input data is in the order of 10^<sup>-16</sup>, often 0.0

The git branch "with-optional-validation-tests" contains the Python and Julia testing infrastructure to compare the correctness of results against the existing FelSiq/DBCV Python implementation, while the main branch includes just a simple CI/CD test with no dependencies on Python. To evaluate the correctness of Dbcv.jl results, checkout the with-optional-validation-tests branch, enter the tests folder, and run ```test_all_datesets.py``` (uses julia/scipy library Kruskal MST) and/or ```test_all_datasets_prim.py``` (uses Dbcv.jl internal PRIM MST implementation and FelSiq/DBCV internal PRIM MST implementation, close to the original MATLAB DBCV code by Pablo A. Jaskowiak: [github.com/pajaskowiak/dbcv](https://github.com/pajaskowiak/dbcv)).

Python modules sklearn (or scikit-learn), numpy, scipy and mpmath are required to be installed in order to be able to run the validation tests in that branch.
One notable difference is in the definition of a custom threshold distance (defaulting to 10^<sup>-9</sup>) below which points are considered duplicates: this threshold is added back to the Validity Validation Score. In Felsiq, this is done adding a fixed quantity (10^<sup>-12</sup>) irrespective of the eventual user set threshold. In Dbcv.jl this values is calculated based on the user set threshold at runtime, which appairs to be more a more correct  behaviour to the author, creating divergences in results with FelSiq/DBCV in the case of a custom set threshold. An optional parameter felsiq_bugforbug can be set to true to reintroduce the broken behaviour in Dbcv.jl for bug for bug compatible results in that case.

The validation suite in python+julia runs on Windows, Mac OS X, GNU/Linux and FreeBSD, and possibly other *nix systems provided the above mentioned python dependecies are available.

Beware that by default FelSiq/DBCV uses Kruskal via Scipy, while we default to using the internal PRIM implementation for closeness with the [original Dbcv paper implementation](https://github.com/pajaskowiak/dbcv). Results between FelSiq/DBCV and Dbcv.jl must be compared between the same MST algorithm. Different Minimum Spanning Tree algorithms, sometimes even different implementations ordering results differently, will largely impact the calculated Dbcv index.

## Installation instructions

Once the packages get, hopefully, registred in Julia Packages, you will be able to install it by running the Julia REPL and invoking
```julia
julia> using Pkg

julia> Pkg.add("DensityBasedClusteringValidation")
```

If not using the Julia package registry or testing this code before its submission request, clone the repository and reference src/DensityBasedClusteringValidation.jl inside your Julia code with its relative or absolute path, like this in the case the src directory is parent of the one containing your own julia sources:

```julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
import DensityBasedClusteringValidation 
```

to simply specify an absolute path irrespective of the current working directory ```@__DIR__```, simply specify the full path without joinpath:

```julia
push!(LOAD_PATH, "path_to_dbcv_src_subfolder")
import DensityBasedClusteringValidation 
```


## Usage instructions

import DensityBasedClusteringValidation in your Julia sources via

```julia
import DensityBasedClusteringValidation 
```

then invoke DensityBasedClusteringValidation.dbcv passing it a dataset and a classification in clusters of that dataset.

```julia
result = DensityBasedClusteringValidation.dbcv(dataset, classification)
```

the two must already be in memory (bound to a variable). The classification vector must only include the classification of the points, in the same order they appear in the dataset, to a cluster id. The dataset must not be repeated in the classification variable. If it is, you can extract a view from your classification matrix, containing only the column with the cluster ids for the original data. If files have to be read first, eg, from csv files, you have to read them into memory before calling DensityBasedClusteringValidation.dbcv, for example via the Julia package DelimitedFiles.

Here follows an example that combines csv reading and extracting only the relevant column (the last one) from the classification file, then chosing to use Kruskal instead of the default PRIM MST:

```julia
import DensityBasedClusteringValidation, DelimitedFiles
has_header::Integer = 0
dataset::AbstractArray = []
dataset_file::String, clustering_file::String = ARGS

if has_header > 0
    dataset = DelimitedFiles.readdlm(dataset_file, ',', BigFloat, skipstart=1)
else
    dataset = DelimitedFiles.readdlm(dataset_file, ',', BigFloat)
end

clustering::AbstractArray = DelimitedFiles.readdlm(clustering_file, ',', Int)

result::Real = DensityBasedClusteringValidation.dbcv(dataset, vec(clustering[:, 1]), use_libgraphs_kruskal=true)

print(result)
```

Note that the code of DensityBasedClusteringValidation.jl is multithreaded, particurarly benefiting from classifications to a large number of clusters (each cluster gets evalued parallely on a separate thread), even if this is not the only implemented parallelism. 
To benefit from multithreading, you must start Julia specifiying the number of available threads to the interpreter/VM, like

```
julia --threads 16
```

or set the environment variable ```JULIA_NUM_THREADS``` on your system.

## Example execution (w/ classification run)

after having installed the DensityBasedClusteringValidation and Clustering Julia packages via pkg: 

```julia
import DensityBasedClusteringValidation, Clustering

points = randn(3, 10000)
clustering = Clustering.dbscan(points, 0.05)

index = DensityBasedClusteringValidation.dbcv(points, clustering.assignments)

print(index)
```


## Supported options

Optional named parameters that can be passed to DensityBasedClusteringValidation.dbcv() are:

```check_duplicates```, defaults to true: whatever to check the dataset for duplicate points (with pairwise distances below the set sep_threshold). Dbcv exits with error if duplicate points are found.

```noise_id```, integer, the cluster id assigned to "noise" points by the classification algorithm used and which is being evaluated. Defaults to considering the '-1' cluster id as the noise cluster id.

```metric```, defaulting to "SqEuclidean" for squared euclidean. The metric is used to calculate pairwise (each point to each point) distances in a cluster and to filter out duplicate points which are closer to each other than the set threshold. Supported metrics are those defined by the Julia package Distances, see the table at https://github.com/JuliaStats/Distances.jl#distance-type-hierarchy
The "type name" as written in Distances.jl docs must be provided as value of the optional named parameter "metric". Note that not all metrics make sense for all datasets and some may just cause errors if incompatible with the data.

```sep_threshold```, the minimum separation distance between points below which to consider points duplicates if check_duplicates is set to true (its default). This value is also indirectly used to clamp (filter) matrices removing smaller values and gets added back, 3 orders of magnitude smaller, to the final result in the step calculating the VCS (Validity Validation Score). Defaults to 1*e^-9.

```felsiq_bugforbug``` FelSiq/DBCV does add back a fixed quantity (1E^-12) in the VCS score, but this does not change when the user sets a custom sep_threshold. We default to deriving that value at runtime based on the user set sep_threshold for correctness, but also implement FelSiq/DBCV behaviour (always add 1E^-12) when the optional named parameter felsiq_bugforbug is set to true. defaults to false.

```bits_of_precision```, defaults to 512 for compatibility with FelSiq/DBCV. It sets the BigFloat type to used a fixed precision instead of dynamically adjusting the precision at runtime, which is otherwise the Julia behaviour. Chosing a bit size compatible with native cpu floating point lenght and instructions will problably speed up execution, reducing the need to split an arithmetic operation over multiple instructions and cpu clock cycle.

```use_libgraphs_kruskal```, default to false. Whatever to use Krustak as the MST to get the minimum spanning tree of points in a cluster, as provided by the Julia Package Libgraphs, or the interal PRIM MST implementation, close to the [original code in MATLAB by the Dbcv paper author](https://github.com/pajaskowiak/dbcv). Note that while FelSiq/DBCV also provides an option in that regarding (internal prim or krustal via scipy), FelSiq/DBCV defaults to using Kruskal.

## Contacts
You can contact the author via issues on [this github repository](https://github.com/kmotavalli/DensityBasedClusteringValidation.jl/issues) or by email at keivan@motavalli.me

This software library has been developed as part of a university internship at the Università di Milano-Bicocca ([www.unimib.it](https://www.unimib.it)) under the supervision of [Davide Chicco](https://www.davidechicco.it) (davide.chicco@unimib.it) of the Department of Informatics, Systems and Communication


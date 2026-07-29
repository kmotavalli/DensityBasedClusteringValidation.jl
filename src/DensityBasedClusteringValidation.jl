module DensityBasedClusteringValidation 

    import Graphs, SimpleWeightedGraphs, Combinatorics, Distances, DataStructures
    using Base.Threads

    export dbcv

    function convert_singleton_clusters_to_noise!(y, noise_id)
        cluster_ids = unique(y)
        cluster_sizes = [count(==(id), y) for id in cluster_ids]
        non_noise_ids::AbstractArray{Bool} = similar(y, Bool)
        non_noise_ids .= true
        for i in eachindex(y)
            pos = findfirst(id -> id == y[i], cluster_ids)
            if cluster_sizes[pos] == 1 || y[i] == noise_id
                #i from the iterator can be a multidimensional coordinates tuple if y is multidim
                y[i] = noise_id
                non_noise_ids[i] = false
            end
        end
        return non_noise_ids
    end

    function dense_ranking(y::AbstractVector)
        vals = sort(unique(y))
        return searchsortedfirst.(Ref(vals), y)
    end

    function subproblem_distances_duplicates!(X::AbstractArray, x_size::Integer, distances::AbstractArray{BigFloat}, dist_instance,
        start_index::Integer, stop_index::Integer, threshold::BigFloat)::Bool
        @inbounds for i in start_index:stop_index
            distances[i, i] = +Inf
            @inbounds for j in i+1:x_size #can also add @simd but not sure if it will propagate to the metric
                dist::BigFloat = dist_instance(view(X, i, :), view(X, j, :))

                if dist < threshold
                    return true
                end

                distances[i, j] = dist
                distances[j, i] = dist
            end
        end
        return false
    end

    function pair_to_pair_distances(X::AbstractArray;
        metric::AbstractString="SqEuclidean",
        threshold::BigFloat=1e-9, dup_check::Bool, felsiq_bugforbug::Bool)::AbstractMatrix

        x_size::Integer = size(X, 1)
        if x_size <= 1
            return nothing
        end
        tolerance::BigFloat = 1e-12
        if !felsiq_bugforbug
            tolerance = threshold / 1000
        end
        distances = Matrix{BigFloat}(undef, x_size, x_size)
        metric_instance = try
        	metric_sym = Symbol(metric)
        	metric_class = getfield(Distances, metric_sym)
		    #instanciate
		    try
			    metric_class(tolerance)
		    catch
			    print("tolerance not supported by " * metric * "\n")
			    metric_class()
		    end
	    catch e
		    error("Metric " * metric * " not found in the Distances.jl package, is it spelled correctly?\n
			    see https://github.com/JuliaStats/Distances.jl?tab=readme-ov-file#distance-type-hierarchy\n")
	    end

        num_threads::Integer = Threads.nthreads(:default)
        problem_size::Integer = div(x_size, num_threads) #integer division
        reminder::Integer = x_size - (problem_size * num_threads)

        if reminder > 0
            num_threads = num_threads + 1
        end

        handles::AbstractVector = Vector{Task}(undef, num_threads)

        cur_index::Integer = 1

        for t in 1:num_threads
            start = cur_index
            if t < num_threads
                finish = cur_index + problem_size
                cur_index = cur_index + problem_size
            else
                if reminder > 0
                    finish = cur_index + reminder -1
                else
                    finish = cur_index + problem_size - 1
                end
            end
            handles[t] = Threads.@spawn subproblem_distances_duplicates!(X, x_size, distances, metric_instance, start, finish, threshold)
        end
        
        results = fetch.(handles)

        for i in 1:num_threads
            if dup_check && results[i] == true
                error("Duplicate samples have been found in X. Try changing sep_threshold (def e^-9)")
                exit()
            end
        end

        clamp!(distances, tolerance, Inf)
        return distances
    end

    function in_cluster_core_distance(distances::AbstractArray{BigFloat},
        d::Integer)::Vector{BigFloat}

        n = size(distances, 1)
        core_dists = (sum(distances .^ -d, dims=ndims(distances)) ./ (n - 1)) .^ (-1.0 / d)

        return vec(core_dists)
    end
    function prim_st(graph::SimpleWeightedGraphs.SimpleWeightedGraph, start::Integer, size::Integer)
        
        cheapestCost::Vector{BigFloat} = fill(Inf, size)
        cheapestEdge::Vector{Union{SimpleWeightedGraphs.SimpleWeightedEdge,Nothing}} = fill(nothing, size)

        explored::Vector{Integer} = []
        bool_explored::Vector{Bool} = fill(false, size)
        unexplored::DataStructures.PriorityQueue = DataStructures.PriorityQueue{Integer, BigFloat}()
        currentVertex::Integer = start

        cheapestCost[start] = 0.0
        DataStructures.enqueue!(unexplored, start, 0.0)

        while !isempty(unexplored)
            currentVertex = DataStructures.dequeue!(unexplored)

            push!(explored, currentVertex)
            bool_explored[currentVertex] = true

            outgoing_edges = [SimpleWeightedGraphs.SimpleWeightedEdge(currentVertex, neighbor, graph.weights[currentVertex, neighbor]) for neighbor in Graphs.outneighbors(graph, currentVertex)]
            for edge in outgoing_edges
                if !bool_explored[edge.dst] && edge.weight < cheapestCost[edge.dst]
                    cheapestCost[edge.dst] = edge.weight 
                    cheapestEdge[edge.dst] = edge

                    unexplored[edge.dst] = edge.weight
                end
            end
        end

        resultEdges::Vector{Tuple{Integer, Integer, BigFloat}} = []
        for vertex in collect(1:size)
            if cheapestEdge[vertex] !== nothing
                edge = cheapestEdge[vertex]
                push!(resultEdges, (edge.src, edge.dst, edge.weight))
            end
        end

        return resultEdges
    end

    function internal_objects(mutual_reach_distances::AbstractMatrix{BigFloat}, use_external_kruskal::Bool)
        
        mrd::AbstractArray{BigFloat} = (mutual_reach_distances + mutual_reach_distances') / 2

        for i in 1:size(mrd, 1)
            mrd[i, i] = 0.0
        end

        graph = SimpleWeightedGraphs.SimpleWeightedGraph(mrd)
        
        n = size(mrd, 1)
        mst_matrix = zeros(eltype(mrd), n, n)

        if use_external_kruskal
            mst_edges = Graphs.kruskal_mst(graph)
            for edge in mst_edges
                src, dest = Graphs.src(edge), Graphs.dst(edge) #these are vertex ids, not floating point coordinates
                w = edge.weight
                mst_matrix[src, dest] = w
                mst_matrix[dest, src] = w
            end
        else
            mst_edges = prim_st(graph, 1, n)
            for edge in mst_edges
                src = edge[1]
                dest = edge[2]
                w = edge[3]
                mst_matrix[src, dest] = w
                mst_matrix[dest, src] = w
            end
        end


        internal_nodes_i = findall(vec(count(>(0.0), mst_matrix, dims=1)) .> 1)
        internal_weights = get_subarray(mst_matrix, internal_nodes_i, nothing)


        if !isempty(internal_nodes_i)
            if length(internal_weights) > 1
                return (internal_nodes_i, internal_weights)
            else
                return (internal_nodes_i, mst_matrix)
            end
        elseif length(internal_weights) > 1
                return (1:size(mutual_reach_distances, 1), internal_weights)
        else
                return (1:size(mutual_reach_distances, 1), mst_matrix)
        end

    end

    function mutual_reachability_distances(mutual_distances::AbstractArray,
        d::Integer)::Tuple{AbstractArray, AbstractArray}

        core_distances::AbstractArray{BigFloat} = in_cluster_core_distance(Matrix(mutual_distances), d)
        cd_appo::AbstractArray{BigFloat} = transpose(core_distances)
        mrd::AbstractArray{BigFloat} = max.(mutual_distances, core_distances, cd_appo)
        return (core_distances, mrd)
    end

    function get_subarray(arr, inds_a, inds_b)
        if isnothing(inds_a)
            return arr
        end

        actual_inds_b = isnothing(inds_b) ? inds_a : inds_b

        return arr[inds_a, actual_inds_b]
    end


    function density_sparseness(cluster_inds::AbstractArray,
        distances::AbstractArray{BigFloat},
        d::Integer, use_external_kruskal::Bool)

        core_distances, mutual_reach_distances = mutual_reachability_distances(distances, d)
        internal_node_inds, internal_edge_weights = internal_objects(mutual_reach_distances, use_external_kruskal)

        valid_weights = internal_edge_weights[isfinite.(internal_edge_weights)]
        cluster_density_sparseness = isempty(valid_weights) ? 0.0 : maximum(valid_weights)
        internal_core_distances = core_distances[internal_node_inds]
        internal_node_inds = cluster_inds[internal_node_inds]

        return (cluster_density_sparseness, internal_core_distances, internal_node_inds)
    end

    function density_separation(cluster_i::Integer,
        cluster_j::Integer,
        distances::AbstractArray{BigFloat},
        internal_core_distances_i::AbstractArray{BigFloat},
        internal_core_distances_j::AbstractArray{BigFloat})


        separation = max.(distances, internal_core_distances_i, transpose(internal_core_distances_j))

        if !isempty(separation)
            density_sep_btwn_clusters = minimum(separation)
        else
            density_sep_btwn_clusters = Inf
        end

        return (cluster_i, cluster_j, density_sep_btwn_clusters)
    end

    function dbcv(Xo::AbstractArray,
        yo::AbstractVector;
        metric::AbstractString = "SqEuclidean",
        noise_id::Integer = -1,
        check_duplicates::Bool = true,
        sep_threshold = 1e-9,
        bits_of_precision = 512,
        use_libgraphs_kruskal::Bool = false,
        felsiq_bugforbug::Bool = false,
    )::AbstractFloat

        setprecision(BigFloat, bits_of_precision)
        if sep_threshold <= 0
            throw(ArgumentError("sep_threshold must be positive"))
        end
        sep_threshold::BigFloat = BigFloat(sep_threshold)

        n, d = size(Xo)

        if n != size(yo, 1)
            Xo = transpose(Xo)
            n, d = size(Xo)
            if n != size(yo, 1)
                throw(ArgumentError("Mismatch in input data (X) lenght and their clustering id assignments (y)"))
            end
        end

        #clusters containing a single element can have that element regarded as noise.
        #reassign all noise to the noise_id cluster, by default id -1

        bool_keep_matrix = convert_singleton_clusters_to_noise!(yo, noise_id)
        
        y = @view yo[bool_keep_matrix]
        #keep whole colums, on true, bool_keep_matrix is not a flattened index for X!
        X = @view Xo[bool_keep_matrix, :]

        if isempty(y)
            return 0.0
        end

        y = dense_ranking(y)

        cluster_ids = unique(y)
        num_clusters::Integer = length(cluster_ids)
        cluster_sizes =  [count(==(id), y) for id in cluster_ids]

        distances::AbstractArray{BigFloat} = pair_to_pair_distances(X, metric=metric, threshold=sep_threshold, dup_check=check_duplicates, felsiq_bugforbug=felsiq_bugforbug)

        # DSC: 'Density Sparseness of a Cluster' init
        dscs = zeros(size(cluster_ids))

        # DSC: 'Density Sparseness of a Cluster' init
        min_dspcs = fill(Inf, size(cluster_ids))

        # core distances of internal nodes
        internal_objects_per_cluster = Vector{AbstractArray}(undef, num_clusters)
        # insert with internal_objects_per_cluster[key] = newvalue

        internal_core_distances_per_cluster = Vector{AbstractArray{BigFloat}}(undef, num_clusters)

        cluster_indexes = [findall(y .== cls_id) for cls_id in cluster_ids]

        Threads.@threads for i in eachindex(cluster_ids)
            subcls_indexes = cluster_indexes[i]
            dscs[i],
            internal_core_distances_per_cluster[i],
            internal_objects_per_cluster[i] =
            density_sparseness(subcls_indexes, get_subarray(distances,subcls_indexes, nothing), d, use_libgraphs_kruskal)
        end

        number_cluster_pairs::Integer = fld((num_clusters*(num_clusters - 1)), 2)

        if number_cluster_pairs > 0
            tasks = [Threads.@spawn density_separation(
                    pair[1],
                    pair[2],
                    get_subarray(distances,
                        internal_objects_per_cluster[pair[1]],
                        internal_objects_per_cluster[pair[2]]),
                    internal_core_distances_per_cluster[pair[1]],
                    internal_core_distances_per_cluster[pair[2]]) 
                for pair in Combinatorics.combinations(cluster_ids, 2)]
            
            results = fetch.(tasks)

            for (cluster_a, cluster_b, density_sep) in results
                min_dspcs[cluster_a] = min(min_dspcs[cluster_a], density_sep)
                min_dspcs[cluster_b] = min(min_dspcs[cluster_b], density_sep)
            end
        end

        # equivalent to np.nan_to_num(min_dspcs, copy=False, posinf=1e12)

        clamp!(min_dspcs, sep_threshold, 1e12)

        base::BigFloat = sep_threshold/1000
        vcs::AbstractArray = []
        if felsiq_bugforbug
            vcs = (min_dspcs .- dscs) ./ (1e-12 .+ max.(min_dspcs, dscs))
        else
            vcs = (min_dspcs .- dscs) ./ (base .+ max.(min_dspcs, dscs))
        end

        # equivalent to np.nan_to_num(vcs, copy=False, nan=0.0)

        replace!(vcs, NaN => 0.0)

        dbcv = sum(vcs .* cluster_sizes) / n
        return dbcv
    end
end

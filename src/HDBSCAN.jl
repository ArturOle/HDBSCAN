

using Distributions
using Plots
using BenchmarkTools
using Base.Threads
using Statistics

include("graphs.jl")


"""
Algorithm 1: HDBSCAN* main steps.
Compute the core distance for the k nearest neighbors for all points in the dataset;
Compute the extended minimum spanning tree from a weighted graph, where the
mutual reachability distances are the edges;
Build the HDBSCAN* hierarchy from the extended minimum spanning tree;
Find the prominent clusters from the hierarchy;
Calculate the outlier scores;
"""
function hdbscan(data, mpts=5)
    (core_distances, precomputed_distances) = precompute_distances(data, mpts)

    #   1) Compute core distance m_pts for all data objects in data
    # core_distances = compute_core_distances(data, precomputed_distances, mpts)
    
    #   2) Compute an MST of G_m_pts - Mutual Reachability Graph
    mrg = build_mutual_reachability_graph(data, precomputed_distances)


    #   3) Extend the MST to obtain MST_ext by adding for each vertex a "self-edge"
    #   with the core distance of corresponding object as weighted

    #   4) Extract the HDBSCAN hierarchy as dendogram from MST_ext

    #   4.1) For the root of the tree assign all objects the same label (single cluster)

    #   4.2) Iteratively remove all edges from MST_ext in decreasing order of weights
    #   (in case of ties, edges must be removed simultaneously)

    #   4.2.1) Before each removal, set the dendogram scale value of the current hierarchical
    #   level as the weight of the edges to be removed

    #   4.2.2) After each removal, assign labels to the connected components that contains the end
    #   vertex of the removed edge to obtain the next hierarchical level: assign a new
    #   cluster label to a component if it still has at least one edge, else assign 
    #   it a null label "noise"

end


function smart_mrg(data, precomp_dist)
    # mrg = graph()
end


function build_mutual_reachability_graph(data, precomp_dist)
    # mrg = Graphs.Graph()
    # display(data)
    # display(precomp_dist)

    # for i in eachindex(data[1])
    # end
end


function precompute_distances(data, m_pts)
    len_data_all_1 = length(data[:, 1])
    distance_matrix = Matrix(undef, len_data_all_1, len_data_all_1)
    core_distances = Vector(undef, len_data_all_1)
    points = Vector{Point}(undef, len_data_all_1)
    paths = Vector{UndirectedPath}(undef, len_data_all_1-1)

    for i in eachindex(data[:, 1])
        
        for j in i:len_data_all_1
            if i==j
                distance_matrix[i, j] = 0
            else
                distance_matrix[i, j] = Euclidian(data[i, :], data[j, :])
                distance_matrix[j, i] = distance_matrix[i, j]
            end

        end
        core_distances[i] = partialsort(distance_matrix[i, :], m_pts+1)
    end

    mrg = Graph(points, path)
    return (core_distances, distance_matrix)
end


function Minkowski(a::Vector, b::Vector, p=2)
    result = 0

    for i in eachindex(a)
        result += abs(a[i]-b[i])^p
    end

    return result^(1/p)
end


function Euclidian(a::Vector, b::Vector)
    result = 0
    
    for i in eachindex(a)
        result += abs2(a[i]-b[i])
    end

    return sqrt(result)
end


function compute_core_distances(data, precomp_dist, m_pts)
    core_distances = Vector(undef, length(data[:, 1]))

    @threads for index in eachindex(data[:, 1])
        core_distances[index] = partialsort(precomp_dist[index, :], m_pts)
    end

    return core_distances
end


function data_gen(nr_of_samples=300)
    data = Matrix{Int64}(undef, nr_of_samples, 2)

    for i=1:floor(Int64, nr_of_samples/3)
        data[i, 1] = floor(Int64, rand(Normal(1 ,40)))
        data[i, 2] = Int(floor(rand(Normal(1 ,40))))

    end
    for i=floor(Int64, nr_of_samples/3):floor(Int64, 2*nr_of_samples/3)
        data[i, 1] = Int(floor(rand(Normal(150 ,20))))
        data[i, 2] = Int(floor(rand(Normal(15 ,50))))

    end
    for i=floor(Int64, 2*nr_of_samples/3):nr_of_samples
        data[i, 1] = Int(floor(rand(Normal(300 ,35))))
        data[i, 2] = Int(floor(rand(Normal(300, 102))))

    end

    return data
end


data = data_gen()
# s = scatter(data[:, 1], data[:, 2])
# display(s)
hdbscan(data)

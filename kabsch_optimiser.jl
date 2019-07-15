# Script that executes the Kabsch algorithm to align two paired sets of points
# by minimising the RMSD between them.
# Takes 2 .xyz files as input
# Program executed as
# 'julia kabsch_optimiser.jl geomA.xyz geomB.xyz'
# Darya Shchepanovska
# 15th July 2019

using DataFrames
using CSV

function read_xyz(input_xyz)
    return CSV.read(input_xyz,ignoreemptylines=true)
end

function translate_to_centroid(coordinate_matrix)
    return coordinate_matrix
end

function cross_covariance_matrix(Pmatrix,Qmatrix)
    return Pmatrix+Qmatrix
end

function optimal_rotation_matrix(CCmatrix)
    return CCmatrix
end

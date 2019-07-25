# Script that executes the Kabsch algorithm to align two paired sets of points
# by minimising the RMSD between them.
# Takes 2 .xyz files as input
# Program executed as
# 'julia kabsch_optimiser.jl geomP.xyz geomQ.xyz'
# Darya Shchepanovska
# 15th July 2019

using DataFrames
using CSV

using LinearAlgebra
using Statistics

function xyz2matrix(input_xyz)
    df = CSV.File(input_xyz,datarow=3,header=false,ignorerepeated=true,delim=' ') |> DataFrame
    matrix = convert(Matrix,df[:,2:4])
    return matrix
end

function translate_to_centroid(coord_matrix)
    center = [mean(coord_matrix[:,1]);mean(coord_matrix[:,2]);mean(coord_matrix[:,3])]
    centroid = transpose(center)
    translated_geom = broadcast(-,coord_matrix,centroid)
    return translated_geom
end

function cross_covariance_matrix(Pmatrix,Qmatrix)
    CCmatrix = transpose(Pmatrix) * Qmatrix
    return CCmatrix
end

function optimal_rotation_matrix(CCmatrix)
    ORmatrix = sqrt(transpose(CCmatrix)*CCmatrix)*inv(CCmatrix)
    return ORmatrix
end

Pgeom = xyz2matrix("test_geometries/hpaldA.xyz")
Qgeom = xyz2matrix("test_geometries/hpaldB.xyz")

normalisedP = (translate_to_centroid(Pgeom))
normalisedQ = (translate_to_centroid(Qgeom))

xcov = cross_covariance_matrix(normalisedP,normalisedQ)

orot = optimal_rotation_matrix(xcov)

display(orot)

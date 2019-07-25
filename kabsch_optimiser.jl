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
    # Imports .xyz file format as data frame, and converts it to an N x 3 Matrix
    df=CSV.File(input_xyz,datarow=3,header=false,ignorerepeated=true,delim=' ')|>DataFrame
    matrix=convert(Matrix,df[:,2:4])
    return matrix
end

function translate_to_centroid(coord_matrix)
    # Normalises the molecular coordinates by centering them.
    center = [mean(coord_matrix[:,1]);mean(coord_matrix[:,2]);mean(coord_matrix[:,3])]
    centroid = transpose(center)
    translated_geom = broadcast(-,coord_matrix,centroid)
    return translated_geom
end

function cross_covariance_matrix(Pmatrix,Qmatrix)
    # Cross covariance matrix gives measure of variability between two matrices
    CCmatrix = transpose(Pmatrix) * Qmatrix
    return CCmatrix
end

function optimal_rotation_matrix(CCmatrix)
    # Returns 3x3 matrix that can be applied to P to get Q
    ORmatrix = sqrt(transpose(CCmatrix)*CCmatrix)*inv(CCmatrix)
    return ORmatrix
end

Pgeom = xyz2matrix(ARGS[1])
Qgeom = xyz2matrix(ARGS[2])

normalisedP = (translate_to_centroid(Pgeom))
normalisedQ = (translate_to_centroid(Qgeom))

xcov = cross_covariance_matrix(normalisedP,normalisedQ)
orot = optimal_rotation_matrix(xcov)

num_atoms = size(normalisedP)[1]
rotated = zeros(Float64,num_atoms,3)

for i=1:num_atoms
    rotated[i,:] = orot*normalisedP[i,:]
end

# display(rotated)

RMSD_value = norm(rotated-normalisedQ)

println("RMSD between these two structures is ", RMSD_value)

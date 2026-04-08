using Gridap, GridapGmsh
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver

using LinearAlgebra, SparseArrays, FillArrays, BlockArrays, WriteVTK

using GridapDistributed, PartitionedArrays
using Logging

using Parameters, TimerOutputs

using MPI, DataFrames, Plots

using DrWatson
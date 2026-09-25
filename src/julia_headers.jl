using Gridap, GridapGmsh
using Gridap.Algebra
using Gridap.CellData
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.MultiField
using Gridap.ReferenceFEs

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers:
  BiformBlock, BlockTriangularSolver, LinearSystemBlock

using BlockArrays
using FillArrays
using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using WriteVTK

using GridapDistributed, PartitionedArrays
using GridapTrilinos
using Logging

using Parameters, TimerOutputs

using MPI

using DrWatson

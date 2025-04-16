module SurfaceBulkViscousFlows

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs

using GridapEmbedded
using GridapEmbedded.AlgoimUtils
using GridapEmbedded.Interfaces: IN, OUT, CUT

using PartitionedArrays
using GridapDistributed

using SparseMatricesCSR
using GridapPETSc
using GridapPETSc: PETSC

using LinearAlgebra: diag, cond
using Symbolics

include("TangentOperators.jl")
include("SolverFunctions.jl")

include("SurfaceBulkInSphere.jl")

export fluid_sphere
export fluid_sphere_axisymmetric
export surface_bulk_in_sphere
export surface_bulk_in_sphere_axisymmetric
export bulk_in_prolate_spheroid
export nuclear_body

end # module SurfaceBulkViscousFlows

module SurfaceBulkViscousFlows

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.FESpaces

using GridapEmbedded
using GridapEmbedded.AlgoimUtils
using GridapEmbedded.Interfaces: IN, OUT, CUT

using SparseMatricesCSR
using GridapPETSc
using GridapPETSc: PETSC

using LinearAlgebra: diag, cond
using Symbolics

import Random

include("ActivityFunctions.jl")
include("AuxiliaryFunctions.jl")
include("InitialDensityFunctions.jl")
include("TangentOperators.jl")
include("SolverFunctions.jl")
include("WeakForms.jl")

include("SurfaceBulkInSphere.jl")
include("SurfaceBulkViscousFlowsAxisymmetric.jl")
include("SurfaceBulkViscousFlows3D.jl")

export unit_density
export verification
export mechanostability

export contractile_ring_axisymmetric
export contractile_ring_3D

export surface_bulk_in_sphere_axisymmetric
export surface_bulk_viscous_flows_axisymmetric
export surface_bulk_viscous_flows_3D

end # module SurfaceBulkViscousFlows

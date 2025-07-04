using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

R = 1.0
n = 25

domain_2D = (-1.2,1.2,0.0,1.2)
domain_3D = (-1.2,1.2,-1.2,1.2,-1.2,1.2)

ls_2D = AlgoimCallLevelSetFunction(
  x -> x[1]^2 + x[2]^4 - 1.0,
  x -> VectorValue( 2.0 * x[1], 4.0 * x[2]^3 ) )

ls_3D = AlgoimCallLevelSetFunction(
  x -> x[1]^2 + x[2]^2 + x[3]^4 - 1.0,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2], 4.0 * x[3]^3 ) )

# Divergence theorem on manifolds (+ algoim quads) allow us to 
# avoid computing the surface divergence of the active stresses 
# (tangent tensor).
#
# In particular, no need to compute mean curvature of the surface.
# See, e.g., (2.18) in Fries and Schöllhammer, "A unified finite 
# strain theory for membranes and ropes", CMAME 365:113031, 2020.
#
# https://arxiv.org/pdf/1909.12640

check_divergence_theorem_on_manifold_2D(
  domain_2D,ls_2D,n,initial_density=verification_2D)

check_divergence_theorem_on_manifold_3D(
  domain_3D,ls_3D,n,initial_density=verification_3D)
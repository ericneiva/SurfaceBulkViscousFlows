using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-1.2,1.2,0.0,1.2)

ls = AlgoimCallLevelSetFunction(
  x -> x[1]*x[1] + x[2]*x[2] - 1.0,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2] ) )

Pe = 30.0
τᵈkₒ = 10.0
μˡ = 1.0e-4
R  = 1.0
n  = 30
Δt = 0.0001
T  = 3.0
output_frequency = 1

GridapPETSc.with() do

  surface_bulk_viscous_flows_axisymmetric(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,initial_density=verification_2D,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name="examples/Verification/plt")

end
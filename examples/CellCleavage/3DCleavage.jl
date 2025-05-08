using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-0.75,0.75,-0.75,0.75,-0.75,0.75)

ls = AlgoimCallLevelSetFunction(
  x -> 4.0 * 1.2 * ( x[1]*x[1] + x[2]*x[2] + x[3]*x[3] ) - 1.0,
  x -> VectorValue( 8.0 * 1.2 * x[1], 8.0 * 1.2 * x[2], 8.0 * 1.2 * x[3] ) )

Pe   = 15.0
τᵈkₒ = 1.0
μˡ   = 1.0e-5
R    = 1.0 / √(4.0*1.2)
n    = 40
Δt   = 0.0001
T    = 1.0
output_frequency = 1

GridapPETSc.with() do

  surface_bulk_viscous_flows_3D(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,γᶜ=10.0,τᵈkₒ=τᵈkₒ,
    initial_density=unit_density,activity=contractile_ring_3D,
    name="examples/CellCleavage/3DCleavage")

end
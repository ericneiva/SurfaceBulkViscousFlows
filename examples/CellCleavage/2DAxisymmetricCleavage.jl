using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-0.75,0.75,0.0,0.75)

ls = AlgoimCallLevelSetFunction(
  x -> 4.0 * 1.2 * ( x[1]*x[1] + x[2]*x[2] ) - 1.0,
  x -> VectorValue( 8.0 * 1.2 * x[1], 8.0 * 1.2 * x[2] ) )

Pe = 5.0
τᵈkₒ = 1.0
μˡ = 1.0e-5
R  = 1.0 / √(4.0*1.2)
n  = 90
Δt = 0.0001
T  = 1.0
output_frequency = 1
redistance_frequency = 1

GridapPETSc.with() do

  surface_bulk_viscous_flows_axisymmetric(
    domain,ls,Pe,μˡ,R,20,Δt,2*Δt,output_frequency=output_frequency,
    writesol=true,γᶜ=10.0,τᵈkₒ=τᵈkₒ,redistance_frequency=redistance_frequency,
    initial_density=unit_density,activity=contractile_ring_axisymmetric,
    name="examples/CellCleavage/2DAxisymmetricCleavage")

  surface_bulk_viscous_flows_axisymmetric(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,γᶜ=10.0,τᵈkₒ=τᵈkₒ,redistance_frequency=redistance_frequency,
    initial_density=unit_density,activity=contractile_ring_axisymmetric,
    name="examples/CellCleavage/2DAxisymmetricCleavage")

end
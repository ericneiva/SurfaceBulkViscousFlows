using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-2.5,2.5,0.0,2.5)

ls = AlgoimCallLevelSetFunction(
  x -> ( x[1]*x[1] + x[2]*x[2] ) - 1.0,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2] ) )

Pe = 50.0
τᵈkₒ = 0.1
μˡ = 1.0e-5
R  = 1.0
n  = 150
Δt = 0.00005
T  = 0.05
output_frequency = 10

GridapPETSc.with() do

  surface_bulk_viscous_flows_axisymmetric(
    domain,ls,Pe,μˡ,R,40,Δt,2*Δt,output_frequency=output_frequency,
    writesol=true,γᶜ=10.0,τᵈkₒ=τᵈkₒ,
    initial_density=initial_ring_axisymmetric,
    activity=unit_activity_axisymmetric,
    name="examples/CellCleavage/2DAxisymmetricCleavage")

  surface_bulk_viscous_flows_axisymmetric(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,γᶜ=10.0,τᵈkₒ=τᵈkₒ,
    initial_density=initial_ring_axisymmetric,
    activity=unit_activity_axisymmetric,
    name="examples/CellCleavage/2DAxisymmetricCleavage")

end
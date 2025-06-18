using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-1.45,1.45,-1.45,1.45)

ls = AlgoimCallLevelSetFunction(
  x -> ( x[1]*x[1] + x[2]*x[2] ) - 1.0,
  x -> VectorValue( x[1], x[2] ) )

Pe = 5.0
τᵈkₒ = 1.0
μˡ = 1.0e-5
R  = 1.0
n  = 80
Δt = 0.0001
T  = 0.01
output_frequency = 1
redistance_frequency = 1

GridapPETSc.with() do

  surface_bulk_viscous_flows_2D(
    domain,ls,Pe,μˡ,R,20,Δt,2*Δt,output_frequency=output_frequency,
    writesol=true,γᶜ=1.0,τᵈkₒ=τᵈkₒ,redistance_frequency=redistance_frequency,
    initial_density=unit_density,activity=contractile_ring_2D,
    name="examples/CellCleavage/2DCleavage")

  surface_bulk_viscous_flows_2D(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,γᶜ=1.0,τᵈkₒ=τᵈkₒ,redistance_frequency=redistance_frequency,
    initial_density=unit_density,activity=contractile_ring_2D,
    name="examples/CellCleavage/2DCleavage")

end
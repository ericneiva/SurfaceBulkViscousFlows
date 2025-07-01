using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-1.22,1.22,-1.22,1.22,-1.22,1.22)

ls = AlgoimCallLevelSetFunction(
  x -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1.0,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2], 2.0 * x[3] ) )

Pe               = 12.5
τᵈkₒ             = 1.0
μˡ               = 1.0e-5
R                = 1.0
n                = 30
Δt               = 0.001
T                = 1.0
center           = VectorValue(0.0,0.0,0.05)
x_axis           = 0.3
y_axis           = 0.3
z_axis           = 0.5
rotation         = 0.0
output_frequency = 10

GridapPETSc.with() do

  surface_bulk_viscous_flows_3D(
    domain,ls,Pe,μˡ,R,20,Δt,2*Δt,center,x_axis,y_axis,z_axis,rotation=rotation,
    output_frequency=output_frequency,writesol=true,
    initial_density=bistable_wave_3D,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name="examples/OocyteAndSpindle/plt")

  surface_bulk_viscous_flows_3D(
    domain,ls,Pe,μˡ,R,n,Δt,T,center,x_axis,y_axis,z_axis,rotation=rotation,
    output_frequency=output_frequency,writesol=true,
    initial_density=bistable_wave_3D,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name="examples/OocyteAndSpindle/plt")

end
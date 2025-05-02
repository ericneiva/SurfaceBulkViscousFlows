using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-1.1,1.1,-1.1,1.1,-1.1,1.1)

function torus(x)
  R = 0.7
  r = 0.3
  x0 = zero(Point{3,typeof(R)})
  _x = x - x0
  (R - sqrt(_x[1]^2+_x[2]^2) )^2 + _x[3]^2 - r^2
end

ls = AlgoimCallLevelSetFunction( x -> torus(x),
                                 x -> ∇(torus,x) )

Pe   = 5.0
τᵈkₒ = 50.0
μˡ   = 1.0e-5
R    = 0.7
n    = 30
Δt   = 0.0001
T    = 2.0
output_frequency = 1

GridapPETSc.with() do

  surface_bulk_viscous_flows_3D(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,initial_density=unit_density,γᶜ=10.0,τᵈkₒ=τᵈkₒ,
    name="examples/RelaxationDynamics/torus")

end
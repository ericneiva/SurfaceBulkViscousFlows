using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-2.25,1.25,-1.75,1.75,-1.75,1.75)

function pear(x)
  R = 2.0
  ( x[1]^2 + √(x[2]*x[2]+x[3]*x[3])^2 ) * 
      ( 1.0 + 2.0*x[1] + 5.0*x[1]^2 + 6.0*x[1]^3 + 6.0*x[1]^4 +
        4.0*x[1]^5 + x[1]^6 - 3.0*√(x[2]*x[2]+x[3]*x[3])^2 - 
        2.0*x[1]*√(x[2]*x[2]+x[3]*x[3])^2 +
        8.0*x[1]^2*√(x[2]*x[2]+x[3]*x[3])^2 + 
        8.0*x[1]^3*√(x[2]*x[2]+x[3]*x[3])^2 + 
        3.0*x[1]^4*√(x[2]*x[2]+x[3]*x[3])^2 + 
        2.0*√(x[2]*x[2]+x[3]*x[3])^4 + 
        4.0*x[1]*√(x[2]*x[2]+x[3]*x[3])^4 + 
        3.0*x[1]^2*√(x[2]*x[2]+x[3]*x[3])^4 + 
        √(x[2]*x[2]+x[3]*x[3])^6 ) - R^2
end

ls = AlgoimCallLevelSetFunction( x -> pear(x),
                                 x -> ∇(pear,x) )

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
    name="examples/RelaxationDynamics/pear")

end
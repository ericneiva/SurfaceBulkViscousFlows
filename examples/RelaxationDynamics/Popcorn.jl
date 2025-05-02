using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-2.01,2.01,-2.01,2.01,-2.01,2.01)

function popcorn(_x)
  r0=0.6
  σ=0.25
  A=2.0
  x0=zero(Point{3,typeof(r0)})
  function point_k(k,r0)
    if 0 <= k && k<=4
      α = 2*k*π/5
      (r0/sqrt(5))*Point(2*cos(α),2*sin(α),1.)
    elseif 5<=k && k <=9
      α = (2*(k-5)-1)*π/5
      (r0/sqrt(5))*Point(2*cos(α),2*sin(α),-1.)
    elseif k==10
      Point(0.,0.,r0)
    else
      Point(0.,0.,-r0)
    end
  end
  x,y,z = _x - x0
  val = sqrt(x^2+y^2+z^2) - r0
  for k in 0:11
    xk,yk,zk = point_k(k,r0)
    val -= A*exp(-((x-xk)^2+(y-yk)^2+(z-zk)^2)/σ^2)
  end
  val
end

ls = AlgoimCallLevelSetFunction( x -> popcorn(x),
                                 x -> ∇(popcorn,x) )

Pe   = 5.0
τᵈkₒ = 50.0
μˡ   = 1.0e-5
R    = 0.6
n    = 35
Δt   = 0.0001
T    = 2.0
output_frequency = 1

GridapPETSc.with() do

  surface_bulk_viscous_flows_3D(
    domain,ls,Pe,μˡ,R,n,Δt,T,output_frequency=output_frequency,
    writesol=true,initial_density=unit_density,γᶜ=10.0,τᵈkₒ=τᵈkₒ,
    name="examples/RelaxationDynamics/popcorn")

end
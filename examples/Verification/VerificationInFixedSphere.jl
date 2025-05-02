using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC
using Plots
using MathTeXEngine
using SurfaceBulkViscousFlows

pythonplot()
default( guidefont  = (14, "serif"),
         tickfont   = (14, "serif"),
         titlefont  = (14, "serif"),
         legendfont = (10, "serif") )

domain = (0.0,2.5,-1.25,1.25)
R = 1.001 # Avoid aligning sphere with mesh edges
ls = AlgoimCallLevelSetFunction(
  x -> ( (x[1]/R)*(x[1]/R) + (x[2]/R)*(x[2]/R) ) - 1.0,
  x -> VectorValue(2.0*(x[1]/(R*R)),2.0*(x[2]/(R*R))) )

γᵅ = 20.0
order = 2

GridapPETSc.with() do

  PETSC.MatMumpsSetIcntl_handle[] == C_NULL && error("MUMPS not available")

  h    = Float64[]
  l2ᵘⁱ = Float64[]; h1ᵘⁱ = Float64[]; l2ᵖⁱ = Float64[]
  l2ᵘˢ = Float64[]; h1ᵘˢ = Float64[]; l2ᵖˢ = Float64[]
  l2ᵘ  = Float64[]; h1ᵘ  = Float64[]; l2ᵖ  = Float64[]

  for n in [20,40,80,160,320,640]

    println("Solving case: ",n)
    _h, _l2ᵘⁱ, _h1ᵘⁱ, _l2ᵖⁱ, _l2ᵘˢ, _h1ᵘˢ, _l2ᵖˢ, _l2ᵘ, _h1ᵘ, _l2ᵖ = 
      surface_bulk_in_sphere_axisymmetric(
        domain,ls,n,order,γᵅ,R,μⁱ=0.1,μˢ=10.0)
    
    @show _l2ᵘⁱ, _h1ᵘⁱ, _l2ᵖⁱ, _l2ᵘˢ, _h1ᵘˢ, _l2ᵖˢ

    push!(h,_h)
    push!(l2ᵘⁱ,_l2ᵘⁱ); push!(h1ᵘⁱ,_h1ᵘⁱ); push!(l2ᵖⁱ,_l2ᵖⁱ)
    push!(l2ᵘˢ,_l2ᵘˢ); push!(h1ᵘˢ,_h1ᵘˢ); push!(l2ᵖˢ,_l2ᵖˢ)
    push!(l2ᵘ ,_l2ᵘ ); push!(h1ᵘ ,_h1ᵘ ); push!(l2ᵖ ,_l2ᵖ )

  end

  plot_name = "surface_bulk_in_sphere_axisymmetric"
  plot(h,[l2ᵘˢ h1ᵘˢ l2ᵖˢ l2ᵘⁱ h1ᵘⁱ l2ᵖⁱ],
      xaxis=:log, yaxis=:log,
      xticks=[1.0,1.0e-1,1.0e-2,1.0e-3],
      xlims=(10^-2.75,10^0.25),
      yticks=[1.0e-2,1.0e-4,1.0e-6,1.0e-8,1.0e-10,1.0e-12],
      ylims=(1.0e-12,1.0e-2),
      xflip=true,
      label=[ L"|| \boldsymbol{U} - \boldsymbol{U}_h {||}_{L^2(\Gamma)}" ;;
              L"|| \boldsymbol{U} - \boldsymbol{U}_h {||}_{H^1(\Gamma)}" ;;
              L"|| P - P_h {||}_{L^2(\Gamma)}" ;;
              L"|| \boldsymbol{u} - \boldsymbol{u}_h {||}_{L^2(\Omega)}" ;;
              L"|| \boldsymbol{u} - \boldsymbol{u}_h {||}_{H^1(\Omega)}" ;;
              L"|| p - p_h {||}_{L^2(\Omega)}" ],
      shape=[:circle :diamond :dtriangle :circle :diamond :dtriangle],
      color=[:blue   :blue    :blue      :orange :orange  :orange],
      style=[:solid  :solid   :solid     :dot    :dot     :dot],
      xlabel="Mesh size",
      ylabel="Error",
      title=string("Surf.-bulk fixed sphere for ",
                   L"\mu_\Gamma = 10.0"," and ",
                   L"\mu_\Omega = 0.1"),
      legend=:bottomleft,
      linewidth=2,
      markersize=6)
  plot!(twinx(),yaxis=:log,ylims=(1.0e-12,1.0e-2),yticks=nothing)
  plot!(twiny(),xaxis=:log,lims=(10^-2.75,10^0.25),xticks=nothing)
  h_sl = [1.0e-1,1.0e-2]
  sl_2 = [5.0e-3,5.0e-5]
  plot!(h_sl,[sl_2],
        style=:dash,
        color=:gray,
        label=L"h^{-2}",
        linewidth=2)
  h_sl = [0.1,1.0e-2]
  sl_3 = [1.0e-8,1.0e-11]
  plot!(h_sl,[sl_3],
        style=:dashdot,
        color=:gray,
        label=L"h^{-3}",
        linewidth=2)
  plot!(size=(600,400),grid=true)
  savefig("examples/Verification/" * plot_name)

end
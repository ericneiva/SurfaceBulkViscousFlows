function check_divergence_theorem_on_manifold_2D(
          domain::Tuple{Vararg{Float64}},
          ls::AlgoimCallLevelSetFunction,
          n::Int;
          initial_density::Function=verification_2D,
          order::Int=2 )

  cells   = (n,div(n,2))
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω       = Triangulation(bgmodel)

  φ₋  = ls
  Vbg = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order))
  _φ₋ = interpolate_everywhere(ls.φ,Vbg)

  degree = order < 3 ? 5 : 2*order # Results depend on integration degree
  squad = Quadrature(algoim,φ₋,degree,phase=CUT)
  s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)

  _,is_nᶜ = narrow_band_triangulation(Ω,_φ₋,Vbg,is_c₋,0.0005)
  Ωᶜ,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_nᶜ,is_c₋)

  N = num_dims(bgmodel)
  reffeʷ = ReferenceFE(lagrangian,VectorValue{N,Float64},order)
  reffeᵉ = ReferenceFE(lagrangian,Float64,order)

  Vʷ = TestFESpace(Ωᶜ,reffeʷ)

  Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
  Uᵉ = TrialFESpace(Vᵉ)

  eₕ = interpolate_everywhere(initial_density(),Uᵉ)

  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  dξ(e) = 4.0 * e / ( 1.0 + e*e ) ^ 2

  nΓ(x) = ls.∇φ(x) / norm(ls.∇φ(x))
  # Doubled mean curvature in axisymmetric coordinates
  twoH(x) = divergence(y->nΓ(y),x) + nΓ(x)⋅iy(x) - nΓ(x)⋅∇(y->nΓ(y),x)⋅(nΓ(x))'

  # writevtk(Ωᶜ,"2D",cellfields=["LS"=>ls.φ,"2H"=>twoH],nsubcells=4)

  f(μ,e) = ∫( ( - ( divᶜ(μ,nΓ) + μ⋅iy ) * (ξ∘(e)) ) * y )dΓ
  g(μ,e) = ∫( ( μ ⋅ ( (dξ∘(e)) * ∇ᵈ(e,nΓ) - twoH * (ξ∘(e)) * nΓ ) ) * y )dΓ

  dμ = get_fe_basis(Vʷ)
  diff = assemble_vector(f(dμ,eₕ),Vʷ) - assemble_vector(g(dμ,eₕ),Vʷ)

  println("2D axisymmetric:")
  println(norm(diff),",",maximum(abs.(diff)),",",minimum(abs.(diff)))

end

function check_divergence_theorem_on_manifold_3D(
          domain::Tuple{Vararg{Float64}},
          ls::AlgoimCallLevelSetFunction,
          n::Int;
          initial_density::Function=verification_3D,
          order::Int=2 )

  cells   = (n,n,n)
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω       = Triangulation(bgmodel)

  φ₋  = ls
  Vbg = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order))
  _φ₋ = interpolate_everywhere(ls.φ,Vbg)

  degree = order < 3 ? 5 : 2*order # Results depend on integration degree
  squad = Quadrature(algoim,φ₋,degree,phase=CUT)
  s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)

  _,is_nᶜ = narrow_band_triangulation(Ω,_φ₋,Vbg,is_c₋,0.0005)
  Ωᶜ,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_nᶜ,is_c₋)

  N = num_dims(bgmodel)
  reffeʷ = ReferenceFE(lagrangian,VectorValue{N,Float64},order)
  reffeᵉ = ReferenceFE(lagrangian,Float64,order)

  Vʷ = TestFESpace(Ωᶜ,reffeʷ)

  Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
  Uᵉ = TrialFESpace(Vᵉ)

  eₕ = interpolate_everywhere(initial_density(),Uᵉ)

  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  dξ(e) = 4.0 * e / ( 1.0 + e*e ) ^ 2

  nΓ(x) = ls.∇φ(x) / norm(ls.∇φ(x))
  # Doubled mean curvature for a 2D surface in 3D
  twoH(x) = divergence(y->nΓ(y),x) # - nΓ(x)⋅∇(y->nΓ(y),x)⋅(nΓ(x))'

  # writevtk(Ωᶜ,"3D",cellfields=["LS"=>ls.φ,"2H"=>twoH],nsubcells=4)

  f(μ,e) = ∫( ( -( divᶜ(μ,nΓ)) * (ξ∘(e)) ) )dΓ
  g(μ,e) = ∫( ( μ ⋅ ( (dξ∘(e)) * ∇ᵈ(e,nΓ) - twoH * (ξ∘(e)) * nΓ ) ) )dΓ

  dμ = get_fe_basis(Vʷ)
  diff = assemble_vector(f(dμ,eₕ),Vʷ) - assemble_vector(g(dμ,eₕ),Vʷ)

  println("3D:")
  println(norm(diff),",",maximum(abs.(diff)),",",minimum(abs.(diff)))

end
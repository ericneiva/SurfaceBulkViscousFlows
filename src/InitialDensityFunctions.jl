unit_density(U,X,Y,dΓ,dΩᶜ,nΓ) = 1.0

function verification(U,X,Y,dΓ,dΩᶜ,nΓ)
  eₐ = 1.0
  Δe = 0.1
  s₀ = 1.0 / 6.0
  x -> ( atan(x[2],-x[1]) > min(s₀,1.0)*pi ) ? eₐ : eₐ + Δe
end

function mechanostability(U,X,Y,dΓ,dΩᶜ,nΓ)
  Random.seed!(1234)
  _eʳ = 0.00001 * randn(Float64,num_free_dofs(U))
  eʳ = FEFunction(U,_eʳ)

  _rᵃ(u,v) = ∫( u*v )dΓ
  _rᵇ(v)   = ∫( eʳ*v )dΓ
  _rᵐ(u,ℓ) = ∫( u*ℓ )dΓ
  _s(u,v) = ∫( 0.1*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  
  # RMK. Stabilisation probably not needed
  rᵃ((u,l),(v,ℓ)) = _rᵃ(u,v) + _rᵐ(u,ℓ) + _rᵐ(v,l) + _s(u,v)
  rᵇ((v,ℓ)) = _rᵇ(v)

  opʳ = AffineFEOperator(rᵃ,rᵇ,X,Y)
  _eₕ,_ = solve(opʳ)
  _eₕ = 1.0 + _eₕ
end

function contraction_wave(U,X,Y,dΓ,dΩᶜ,nΓ)
  eₐ = 0.8; Δe = 0.4; s₀ = 0.0;
  (x,s) -> ( atan(x[2],x[1]) > min(max(0.0,(s+s₀)*(s+s₀)),0.9*0.9)*pi ) ? eₐ + Δe : eₐ
end
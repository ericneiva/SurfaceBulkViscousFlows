unit_density(U,X,Y,dΓ,dΩᶜ,nΓ) = 1.0

function verification(U,X,Y,dΓ,dΩᶜ,nΓ)
  eₐ = 1.0
  Δe = 0.1
  s₀ = 1.0 / 6.0
  x -> ( atan(x[2],-x[1]) > min(s₀,1.0)*pi ) ? eₐ : eₐ + Δe
end

function mechanostability_axisymmetric(U,X,Y,dΓ,dΩᶜ,nΓ)
  Random.seed!(1234)
  _eʳ = 0.00001 * randn(Float64,num_free_dofs(U))
  eʳ = FEFunction(U,_eʳ)

  y(x) = x[2]

  _rᵃ(u,v) = ∫( (u*v)*y )dΓ
  _rᵇ(v)   = ∫( (eʳ*v)*y )dΓ
  _rᵐ(u,ℓ) = ∫( (u*ℓ)*y )dΓ
  _s(u,v) = ∫( 0.1*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  
  # RMK. Stabilisation probably not needed
  rᵃ((u,l),(v,ℓ)) = _rᵃ(u,v) + _rᵐ(u,ℓ) + _rᵐ(v,l) + _s(u,v)
  rᵇ((v,ℓ)) = _rᵇ(v)

  opʳ = AffineFEOperator(rᵃ,rᵇ,X,Y)
  _eₕ,_ = solve(opʳ)
  _eₕ = 1.0 + _eₕ
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

function initial_ring_axisymmetric(U,X,Y,dΓ,dΩᶜ,nΓ)
  Am = 1.0; Aₒ = 0.05; ω = 0.1;
  x -> ( Aₒ + ( Am - Aₒ ) * ( exp( -x[1]*x[1] / ( 2.0 * ω * ω ) ) ) ) * x[2]
end

function initial_ring_3D(U,X,Y,dΓ,dΩᶜ,nΓ)
  Am = 10.0; Aₒ = 1.0; ω = 0.1;
  x -> ( Aₒ + ( Am - Aₒ ) * ( exp( -x[1]*x[1] / ( 2.0 * ω * ω ) ) ) )
end
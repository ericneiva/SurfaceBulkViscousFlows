iy(x) = VectorValue( 0.0, 1.0 / x[2] ); y(x) = x[2]

function cortical_flow_problem_axisymmetric(ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,
    γ::Float64,Pe::Float64,μˡ::Float64,R::Float64,ξ₀::Function)

  aʷ(υ,μ) = ∫( 2.0 * ( εᶜ(υ,nΓ)⊙εᵈ(μ,nΓ) + (υ⋅iy)*(μ⋅iy) ) * y )dΓ

  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )

  f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ)+μ⋅iy)*(ξ∘(e)) ) * ξ₀ )dΓ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  βʳ(μ,u,p) = ∫( ( μ⋅((σᵘ∘(ε(u),p))⋅nΓ) ) * y )dΓ

  sᵘ(υ,μ) = ∫( γ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  RB¹ = VectorValue(1.0,0.0)

  r¹(u,ℓ) = ∫( ( u⋅(ℓ*RB¹) )*y )dΓ
  r²(u,ℓ) = ∫( ( u⋅(ℓ*nΓ ) )*y )dΓ

  aᵛ((υ,l¹,l²),(μ,ℓ¹,ℓ²)) =
    aʷ(υ,μ) + sᵘ(υ,μ) + r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²)
  bᵛ((μ,ℓ¹,ℓ²)) = f(μ,eₕ) - βʳ(μ,ulₕ,plₕ)

  aᵛ, bᵛ
end

function bulk_flow_problem_axisymmetric(
    υ,dΩ,dΓ,nΓ,μˡ::Float64,R::Float64,γ::Float64,h::Float64)
  
  aᵇ(u,v) = ∫( 2.0 * μˡ * R * (ε(u)⊙ε(v) + (u⋅iy)*(v⋅iy)) * y )dΩ

  bᵇ(v,q) = ∫( q*(∇⋅v + v⋅iy)*y )dΩ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  αˡ(u,v,p,q) = ∫( ( (γ/h)*(u⋅v)           -
                      u⋅((σᵘ∘(ε(v),q))⋅nΓ) -
                      v⋅((σᵘ∘(ε(u),p))⋅nΓ) )*y )dΓ
  αʳ(v,q,υ)   = ∫( ( (γ/h)*(υ⋅v) - υ⋅((σᵘ∘(ε(v),q))⋅nΓ) )*y )dΓ

  r(p,ℓ) = ∫( ( p*ℓ ) * y )dΩ

  aᵘ((uˡ,pˡ,l),(vˡ,qˡ,ℓ)) = 
    aᵇ(uˡ,vˡ) - bᵇ(vˡ,pˡ) - bᵇ(uˡ,qˡ) + 
    αˡ(uˡ,vˡ,pˡ,qˡ) + r(pˡ,ℓ) + r(qˡ,l)
  bᵘ((vˡ,qˡ,l)) = αʳ(vˡ,qˡ,υ)

  aᵘ, bᵘ
end

function transport_problem_axisymmetric(u,eₕ,dΓ,dΩᶜ,nΓ,
    dt::Float64,γ::Float64,τᵈkₒ::Float64)
  
  m(e,ε)  = ∫( (1/dt)*(e*ε)*y )dΓ
  sᵈ(e,ε) = ∫( ( ∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ) )*y )dΓ
  c(e,ε)  = ∫( ( (u⋅∇ᵈ(e,nΓ))*ε + (tr(∇ᵈ(u,nΓ))+u⋅iy)*(e*ε) )*y )dΓ
  r(e,ε)  = ∫( τᵈkₒ*(e*ε)*y )dΓ
  l(ε)    = ∫( τᵈkₒ*ε*y )dΓ
  s(υ,μ)  = ∫( γ*((nΓ⋅∇(υ))⊙(nΓ⋅∇(μ))) )dΩᶜ

  aᵉ(e,ε) = m(e,ε) + c(e,ε) + r(e,ε) + sᵈ(e,ε) + s(e,ε)
  bᵉ(ε)   = m(eₕ,ε) + l(ε)

  aᵉ,bᵉ
end

function transport_problem_axisymmetric(u,eₕ,dΓ,nΓ,γⁿ::Float64,
    dΩᶜ,γ¹::Float64,dΛ,nΛ,τᵈkₒ::Float64,dt::Float64)
  
  m(e,ε) = ∫( (1/dt)*(e*ε)*y )dΓ
  d(e,ε) = ∫( ( ∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ) )*y )dΓ
  c(e,ε) = ∫( ( (u⋅∇ᵈ(e,nΓ))*ε + (tr(∇ᵈ(u,nΓ))+u⋅iy)*(e*ε) )*y )dΓ
  r(e,ε) = ∫( τᵈkₒ*(e*ε)*y )dΓ
  l(ε)   = ∫( τᵈkₒ*ε*y )dΓ
  
  sⁿ(υ,μ) = ∫( γⁿ*((nΓ⋅∇(υ))⊙(nΓ⋅∇(μ))) )dΩᶜ
  s¹(υ,μ) = ∫( γ¹*(jump(nΛ⋅∇(υ))*jump(nΛ⋅∇(μ))) )dΛ

  aᵉ(e,ε) = m(e,ε) + c(e,ε) + r(e,ε) + d(e,ε) + sⁿ(e,ε) + s¹(e,ε)
  bᵉ(ε)   = m(eₕ,ε) + l(ε)

  aᵉ,bᵉ
end

function cortical_flow_problem_3D(ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,
    γ::Float64,Pe::Float64,μˡ::Float64,R::Float64,ξ₀::Function)

  β = 0.001
  aʷ(υ,μ) = ∫( 2.0 * ( εᶜ(υ,nΓ)⊙εᵈ(μ,nΓ) ) + β * (υ⋅μ) )dΓ

  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )

  f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ))*(ξ∘(e)) ) * ξ₀ )dΓ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  βʳ(μ,u,p) = ∫( μ⋅((σᵘ∘(ε(u),p))⋅nΓ) )dΓ

  sᵘ(υ,μ) = ∫( γ*((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  RB¹ = VectorValue(1.0,0.0,0.0)
  RB² = VectorValue(0.0,1.0,0.0)
  RB³ = VectorValue(0.0,0.0,1.0)
  RB⁴ = x -> VectorValue(0.0,-x[3],x[2])
  RB⁵ = x -> VectorValue(x[3],0.0,-x[1])
  RB⁶ = x -> VectorValue(-x[2],x[1],0.0)

  r¹(u,ℓ) = ∫( (u⋅RB¹)*ℓ )dΓ
  r²(u,ℓ) = ∫( (u⋅RB²)*ℓ )dΓ
  r³(u,ℓ) = ∫( (u⋅RB³)*ℓ )dΓ
  r⁴(u,ℓ) = ∫( (u⋅RB⁴)*ℓ )dΓ
  r⁵(u,ℓ) = ∫( (u⋅RB⁵)*ℓ )dΓ
  r⁶(u,ℓ) = ∫( (u⋅RB⁶)*ℓ )dΓ
  r⁷(u,ℓ) = ∫( (u⋅nΓ )*ℓ )dΓ

  aᵛ((υ,l¹,l²,l³,l⁴,l⁵,l⁶,l⁷),(μ,ℓ¹,ℓ²,ℓ³,ℓ⁴,ℓ⁵,ℓ⁶,ℓ⁷)) =
    aʷ(υ,μ) + sᵘ(υ,μ) + 
    r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²) +
    r³(υ,ℓ³) + r³(μ,l³) + r⁴(υ,ℓ⁴) + r⁴(μ,l⁴) +
    r⁵(υ,ℓ⁵) + r⁵(μ,l⁵) + r⁶(υ,ℓ⁶) + r⁶(μ,l⁶) +
    r⁷(υ,ℓ⁷) + r⁷(μ,l⁷)
  bᵛ((μ,ℓ¹,ℓ²,ℓ³,ℓ⁴,ℓ⁵,ℓ⁶,ℓ⁷)) = f(μ,eₕ) - βʳ(μ,ulₕ,plₕ)

  aᵛ, bᵛ
end

function bulk_flow_problem_3D(
  υ,dΩ,dΓ,nΓ,μˡ::Float64,R::Float64,γ::Float64,h::Float64)

  aᵇ(u,v) = ∫( 2.0 * μˡ * R * (ε(u)⊙ε(v)) )dΩ

  bᵇ(v,q) = ∫( q*(∇⋅v) )dΩ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  αˡ(u,v,p,q) = ∫( (γ/h)*(u⋅v)           -
                    u⋅((σᵘ∘(ε(v),q))⋅nΓ) -
                    v⋅((σᵘ∘(ε(u),p))⋅nΓ) )dΓ
  αʳ(v,q,υ) = ∫( (γ/h)*(υ⋅v) - υ⋅((σᵘ∘(ε(v),q))⋅nΓ) )dΓ

  r(p,ℓ) = ∫( p*ℓ )dΩ

  aᵘ((uˡ,pˡ,l),(vˡ,qˡ,ℓ)) = 
    aᵇ(uˡ,vˡ) - bᵇ(vˡ,pˡ) - bᵇ(uˡ,qˡ) + 
    αˡ(uˡ,vˡ,pˡ,qˡ) + r(pˡ,ℓ) + r(qˡ,l)
  bᵘ((vˡ,qˡ,l)) = αʳ(vˡ,qˡ,υ)

  aᵘ, bᵘ
end

function transport_problem_3D(u,eₕ,dΓ,dΩᶜ,nΓ,
    dt::Float64,γ::Float64,τᵈkₒ::Float64)

  m(e,ε)  = ∫( (1/dt)*(e*ε) )dΓ
  sᵈ(e,ε) = ∫( ∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ) )dΓ
  c(e,ε)  = ∫( (u⋅∇ᵈ(e,nΓ))*ε + tr(∇ᵈ(u,nΓ))*(e*ε) )dΓ
  r(e,ε)  = ∫( τᵈkₒ*(e*ε) )dΓ
  l(ε)    = ∫( τᵈkₒ*ε )dΓ
  s(υ,μ)  = ∫( γ*((nΓ⋅∇(υ))⊙(nΓ⋅∇(μ))) )dΩᶜ

  aᵉ(e,ε) = m(e,ε) + c(e,ε) + r(e,ε) + sᵈ(e,ε) + s(e,ε)
  bᵉ(ε)   = m(eₕ,ε) + l(ε)

  aᵉ,bᵉ
end
function bulk_flow_axisymmetric(U::Float64,a::Float64,μ::Float64)

  r(x) = √( x[1]*x[1] + x[2]*x[2] )
  s(x) = x[1] / r(x)
  c(x) = x[2] / r(x)

  eʳ(x) = VectorValue( s(x),  c(x) )
  eᶿ(x) = VectorValue( c(x), -s(x) )

  σ = 1.0
  E =   (5/2) * ( U / ( a * a ) ) * σ / ( 1 + σ )
  G = - (1/4) *   U               * σ / ( 1 + σ )

  uˣ(x) = x[1] * x[2] * E/5                                 # uʳ * s + uᶿ * c
  uᶻ(x) = - (E/5) * ( 2.0*x[1]*x[1] + x[2]*x[2] ) - 2.0 * G # uʳ * c - uᶿ * s

  u(x) = VectorValue(uˣ(x),uᶻ(x))
  f(x) = VectorValue(0.0,2.0*μ*E)
  g(x) = 0.0

  ε(x) = 0.5 * ( ∇(u)(x) + transpose(∇(u)(x)) )
  n(x) = VectorValue( x[1] / r(x), x[2] / r(x) )
  tˢ(x) = 2.0 * μ * (ε(x)⋅n(x))

  u,f,g,tˢ

end

function surface_flow_axisymmetric(U::Float64,a::Float64,μ::Float64)

  @variables x,z

  n = [x/a, z/a]

  σ = 1.0
  E =   (5/2) * ( U / ( a * a ) ) * σ / ( 1 + σ )
  G = - (1/4) *   U               * σ / ( 1 + σ )

  uˣ = x * z * E/5
  uᶻ = - (E/5) * ( 2.0*x*x + z*z ) - 2.0 * G

  _u = [uˣ,uᶻ]
  
  n_tensor_n = n*transpose(n)
  P = one(n_tensor_n) - n_tensor_n
  
  ∇u = Symbolics.jacobian(_u,[x,z])
  
  ∇ᶜu = P*∇u*P
  _g = sum(diag(∇ᶜu)) + uˣ/x

  ∂ˣuˣ = ∇ᶜu[1,1]
  ∂ˣuᶻ = ∇ᶜu[2,1]

  εᶜu = 0.5*(∇ᶜu+transpose(∇ᶜu))

  ᵃεᶜuˣ = 
    sum(diag(P*Symbolics.jacobian(εᶜu[1,:],[x,z]))) + ∂ˣuˣ/x - uˣ/(x^2)
  ᵃεᶜuᶻ = 
    sum(diag(P*Symbolics.jacobian(εᶜu[2,:],[x,z]))) + ∂ˣuᶻ/x
  
  _f₁ = -2.0 * μ * ᵃεᶜuˣ
  _f₂ = -2.0 * μ * ᵃεᶜuᶻ
  _f = P * [_f₁,_f₂]

  u_s = build_function(_u,[x,z],expression=Val{false})[1]
  u(q) = VectorValue(u_s(get_array(q)))

  f_s = build_function(_f,[x,z],expression=Val{false})[1]
  f(q) = VectorValue(f_s(get_array(q)))

  g_s = build_function(_g,[x,z],expression=Val{false})
  g(q) = g_s(get_array(q))

  u,f,g

end

function surface_bulk_in_sphere_axisymmetric(domain::NTuple{D,Float64},
                                             φ::AlgoimCallLevelSetFunction,
                                             n::Int,
                                             order::Int,
                                             γᵅ::Float64,
                                             R::Float64;
                                             partition::Tuple=(n,n),
                                             μⁱ::Float64=1.0,
                                             μˢ::Float64=1.0,
                                             writesol::Bool=true) where {D}

  # Background geometry
  h = (domain[2]-domain[1])/n
  bgmodel = CartesianDiscreteModel(domain,partition)
  Ω = Triangulation(bgmodel)

  # Active triangulations, measures and normal
  degree = order < 3 ? 3 : 2*order
  squad = Quadrature(algoim,φ,degree,phase=CUT)
  Ωᶜ,dΓ,is_c = TriangulationAndMeasure(Ω,squad)
  dΩᶜ = Measure(Ωᶜ,2*order)
  vquad = Quadrature(algoim,φ,degree,phase=IN)
  Ωⁱ,dΩⁱ,is_i = TriangulationAndMeasure(Ω,vquad)
  nΓ = normal(φ,Ω)

  # Reference FEs
  N = num_dims(bgmodel)
  reffeᵘ = ReferenceFE(lagrangian,VectorValue{N,Float64},order)
  reffeᵉ = ReferenceFE(lagrangian,VectorValue{N,Float64},order,space=:S)
  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1,space=:P)
  reffeᵗ = ReferenceFE(lagrangian,Float64,order-1)

  # FE spaces
  cell_to_cellin  = aggregate(Ω,is_i,is_c,IN)

  Vⁱstdᵘ = TestFESpace(Ωⁱ,reffeᵘ,dirichlet_tags=[7])
  Vⁱserᵘ = TestFESpace(Ωⁱ,reffeᵉ,conformity=:L2)
  Vⁱstdᵖ = TestFESpace(Ωⁱ,reffeᵖ)

  Vⁱ = AgFEMSpace(Vⁱstdᵘ,cell_to_cellin,Vⁱserᵘ)
  Qⁱ = AgFEMSpace(Vⁱstdᵖ,cell_to_cellin)

  Vˢ = TestFESpace(Ωᶜ,reffeᵘ,dirichlet_tags=[7])
  Qˢ = TestFESpace(Ωᶜ,reffeᵗ)

  Tˡ = ConstantFESpace(bgmodel)

  U = 1.0
  uⁱ,fⁱ,gⁱ,tˢ = bulk_flow_axisymmetric(U,R,μⁱ)
  uˢ,fˢ,gˢ = surface_flow_axisymmetric(U,R,μˢ)

  Uⁱ = TrialFESpace(Vⁱ,uⁱ)
  Pⁱ = TrialFESpace(Qⁱ)

  Uˢ = TrialFESpace(Vˢ,uˢ)
  Pˢ = TrialFESpace(Qˢ)

  Sˡ = TrialFESpace(Tˡ)

  Yᵘ = MultiFieldFESpace([Vⁱ,Qⁱ,Tˡ])
  Xᵘ = MultiFieldFESpace([Uⁱ,Pⁱ,Sˡ])

  Yˢ = MultiFieldFESpace([Vˢ,Qˢ,Tˡ])
  Xˢ = MultiFieldFESpace([Uˢ,Pˢ,Sˡ])

  # *** BULK FLOW ***
  r(x) = VectorValue( 1.0/x[1], 0.0 )
  dr(x) = x[1]
  # ** grad-grad term **
  aᵘ(u,v,dΩ) = ∫( 2.0*μⁱ * ( ε(u)⊙ε(v) + (u⋅r)*(v⋅r) ) * dr )dΩ
  # ** pressure term **
  bᵘ(v,q,dΩ) = ∫( ( q*(∇⋅v) + q*(v⋅r) ) * dr )dΩ
  # ** Nitsche & coupling terms **
  σᵘ(ε,q) = 2.0*μⁱ*ε - q*one(ε)
  αˡ(u,v,p,q) = ∫( ( (γᵅ/h)*(u⋅v)          -
                      u⋅((σᵘ∘(ε(v),q))⋅nΓ) -
                      v⋅((σᵘ∘(ε(u),p))⋅nΓ) ) * dr )dΓ
  αʳ(v,q,υ) = ∫( ( (γᵅ/h)*(υ⋅v) - υ⋅((σᵘ∘(ε(v),q))⋅nΓ) ) * dr )dΓ
  # ** mean-value constraint **
  sᵖ(p,ℓ,dΩ) = ∫( ( p*ℓ ) * dr )dΩ
  # ** source terms **
  f(v,f,dΩ) = ∫( ( f⋅v ) * dr )dΩ
  g(q,g,dΩ) = ∫( ( g⋅q ) * dr )dΩ

  # *** SURFACE FLOW ***
  # ** mass term **
  m(u,v) = ∫( ( xₜ(u,nΓ)⋅xₜ(v,nΓ) ) * dr )dΓ
  # ** grad-grad term **
  _nΓ(x) = φ.∇φ(x)/norm(φ.∇φ(x))
  H(x) = ∇(_nΓ,x) # TO-CLEAN
  aᵛ(u,v,μ) = ∫( 2.0 * μ * ( (xₜ(u,nΓ)⋅r)*(xₜ(v,nΓ)⋅r) +
    ((εᶜ(u,nΓ)-(u⋅nΓ)*H)⊙(εᶜ(v,nΓ)-(v⋅nΓ)*H)) ) * dr )dΓ
  # ** pressure term **
  bᵛ(u,q) = ∫( ( u⋅∇ᵈ(q,nΓ) ) * dr )dΓ
  # ** weak tangentiality **
  η = 10.0 * μˢ / h^order # 1.0 / (h^(order)) also converges (right quantities?)
  k(u,v) = ∫( η * ((u⋅nΓ)*(v⋅nΓ)) )dΓ
  # ** u- and p-stabilisation **
  γʷ = 0.1 * μˢ / h
  γᵖ = 0.1 * h
  s(u,v,γ) = ∫( γ*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  # ** bulk to surface load **
  βʳ(μ,u,p) = ∫( ( μ⋅((σᵘ∘(ε(u),p))⋅nΓ) ) * dr )dΓ

  uⁱₕ = zero(Uⁱ)
  pⁱₕ = zero(Pⁱ)
  uˢₕ = zero(Uˢ)
  pˢₕ = zero(Pˢ)

  aᵘ((uⁱ,pⁱ,lⁱ),(vⁱ,qⁱ,ℓⁱ)) = 
    aᵘ(uⁱ,vⁱ,dΩⁱ) - bᵘ(vⁱ,pⁱ,dΩⁱ) - bᵘ(uⁱ,qⁱ,dΩⁱ) + 
    sᵖ(pⁱ,ℓⁱ,dΩⁱ) + sᵖ(qⁱ,lⁱ,dΩⁱ) + αˡ(uⁱ,vⁱ,pⁱ,qⁱ)
  bᵘ((vⁱ,qⁱ,ℓⁱ)) =  f(vⁱ,fⁱ,dΩⁱ) - g(qⁱ,gⁱ,dΩⁱ) + αʳ(vⁱ,qⁱ,uⁱ)

  aˢ((uˢ,pˢ,lˢ),(vˢ,qˢ,ℓˢ)) = 
    aᵛ(uˢ,vˢ,μˢ) + bᵛ(vˢ,pˢ) + bᵛ(uˢ,qˢ) + 
    s(uˢ,vˢ,γʷ) - s(pˢ,qˢ,γᵖ) + k(uˢ,vˢ) +
    sᵖ(pˢ,ℓˢ,dΓ) + sᵖ(qˢ,lˢ,dΓ)
  bˢ((vˢ,qˢ,ℓˢ)) = f(vˢ,fˢ,dΓ) - g(qˢ,gˢ,dΓ) + 
                   f(vˢ,tˢ,dΓ) - βʳ(vˢ,uⁱₕ,pⁱₕ)

  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  ps = PETScLinearSolver(mykspsetup)
  assemᵘ = SparseMatrixAssembler(Tm,Tv,Xᵘ,Yᵘ)
  assemˢ = SparseMatrixAssembler(Tm,Tv,Xˢ,Yˢ)

  maxiter = 10
  reltol = 1e-7
  for i in 1:maxiter

    _uˢₕ = uˢₕ
    _uⁱₕ = uⁱₕ
    opˢ = AffineFEOperator(aˢ,bˢ,Xˢ,Yˢ,assemˢ)
    uˢₕ,pˢₕ = solve_stokes(opˢ,Xˢ,ps)
    opᵘ = AffineFEOperator(aᵘ,bᵘ,Xᵘ,Yᵘ,assemᵘ)
    uⁱₕ,pⁱₕ = solve_stokes(opᵘ,Xᵘ,ps)

    chk1 = √( ∑( ∫( (uˢₕ-_uˢₕ)⋅(uˢₕ-_uˢₕ) )dΓ ) ) / √( ∑( ∫( _uˢₕ⋅_uˢₕ )dΓ ) )
    chk2 = √( ∑( ∫( (uⁱₕ-_uⁱₕ)⋅(uⁱₕ-_uⁱₕ) )dΓ ) ) / √( ∑( ∫( _uⁱₕ⋅_uⁱₕ )dΓ ) )
    chk = max(chk1,chk2)

    @info "chk1 = $chk1 and chk2 = $chk2 at i = $i"
    if chk < reltol
      break
    elseif i == maxiter
      error("No convergence")
    end

  end

  l2ᵇ(w,dΩ) = √( ∑( ∫( ( w⊙w ) * dr )dΩ ) )
  h1ᵇ(w,dΩ) = √( ∑( ∫( ( ε(w)⊙ε(w) + (w⋅r)*(w⋅r) ) * dr )dΩ ) )

  eⁱ = uⁱ - uⁱₕ
  il2ᵘ = l2ᵇ(eⁱ,dΩⁱ)
  ih1ᵘ = h1ᵇ(eⁱ,dΩⁱ)
  il2ᵖ = l2ᵇ(pⁱₕ,dΩⁱ)

  l2ˢ(w,dΓ) = √( ∑( ∫( ( w⊙w ) * dr )dΓ ) )
  h1ˢ(w,dΓ,nΓ) = √( ∑( ∫( ( εᶜ(w,nΓ)⊙εᶜ(w,nΓ) + (w⋅r)*(w⋅r) ) * dr )dΓ ) )

  eˢ = uˢ - uˢₕ
  sl2ᵘ = l2ˢ(eˢ,dΓ)
  sh1ᵘ = h1ˢ(eˢ,dΓ,nΓ)
  sl2ᵖ = l2ˢ(pˢₕ,dΓ)

  l2ᵘ = √( il2ᵘ*il2ᵘ + sl2ᵘ*sl2ᵘ )
  h1ᵘ = √( ih1ᵘ*ih1ᵘ + sh1ᵘ*sh1ᵘ )
  l2ᵖ = √( il2ᵖ*il2ᵖ + sl2ᵖ*sl2ᵖ )

  h, sl2ᵘ, sh1ᵘ, sl2ᵖ, il2ᵘ, ih1ᵘ, il2ᵖ, l2ᵘ, h1ᵘ, l2ᵖ

end
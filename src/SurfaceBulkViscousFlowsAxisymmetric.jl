function surface_bulk_viscous_flows_axisymmetric(
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            μˡ::Float64,
            R::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64;
            initial_density::Function = verification,
            activity::Function = unit_activity_axisymmetric,
            order::Int = 2,
            γᶜ::Float64 = 1.0,
            γᵅ::Float64 = 20.0,
            τᵈkₒ::Float64 = 10.0,
            writesol::Bool = true,
            reltol::Float64 = 1.0e-5,
            maxiter::Int = 10,
            output_frequency::Int = 1,
            redistance_frequency::Int = 1,
            name::String = "plt")

  # Background geometry
  cells   = (n,div(n,2))
  h       = (domain[2]-domain[1])/n
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω       = Triangulation(bgmodel)

  # Buffer of active model and integration objects
  degree = order < 3 ? 3 : 2*order
  buffer = Ref{Any}(( Ωᶜ  = nothing, dΩᶜ   = nothing,
                      Ωˡ  = nothing, dΩˡ   = nothing,
                      dΓ  = nothing, nΓ    = nothing,
                      φ₋  = nothing, aggsˡ = nothing,
                      cp₋ = nothing, t     = nothing,
                      Vbg = nothing ))

  function update_buffer!(i,t,dt,v₋₂,mv₋₂)

    if buffer[].t == t
      return true
    else

      Ωᶜ    = buffer[].Ωᶜ
      dΩᶜ   = buffer[].dΩᶜ
      Ωˡ    = buffer[].Ωˡ
      dΩˡ   = buffer[].dΩˡ
      aggsˡ = buffer[].aggsˡ
      dΓ    = buffer[].dΓ
      nΓ    = buffer[].nΓ
      cp₋   = buffer[].cp₋
      φ₋    = buffer[].φ₋
      t     = buffer[].t
      Vbg   = buffer[].Vbg

      if buffer[].Ωᶜ === nothing
        Vbg = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order))
        _φ₋ = interpolate_everywhere(ls.φ,Vbg)
      else
        cp₋₂ = buffer[].cp₋
        φ₋₂  = buffer[].φ₋
        __φ  = get_free_dof_values(φ₋₂.φ)
        Ωⱽ   = get_triangulation(Vbg)
        _ϕ₋  = compute_normal_displacement(cp₋₂,φ₋₂,v₋₂,dt,Ωⱽ)
        ϕ₋   = __φ - _ϕ₋
        _φ₋  = FEFunction(Vbg,ϕ₋)
      end

      # Current time level set
      φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      ( i % redistance_frequency == 0 ) && begin
        _φ₋ = compute_distance_fe_function(bgmodel,Vbg,φ₋,order,cppdegree=3)
        φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      end

      cp₋ = compute_closest_point_projections(Vbg,φ₋,order,
              cppdegree=3,trim=true,limitstol=1.0e-2)

      # Current time surface and bulk measures
      squad = Quadrature(algoim,φ₋,degree,phase=CUT)
      s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)
      vquad = Quadrature(algoim,φ₋,degree,phase=IN)
      v_cell_quad,is_a₋ = CellQuadratureAndActiveMask(bgmodel,vquad)

      # Surface narrow-band triangulation
      δ₋ = 2.0 * mv₋₂ * dt
      _,is_nᶜ = narrow_band_triangulation(Ω,_φ₋,Vbg,is_c₋,δ₋)

      # Narrow band bulk active triangulations and measures
      # For bulk, only need to extend along positive LS vals
      _,is_nᵃ = active_triangulation(Ω,_φ₋,Vbg,is_a₋,δ₋)

      # Current aggregates
      aggsˡ = aggregate(Ω,is_a₋,is_c₋,IN,is_nᵃ)

      Ωˡ,dΩˡ = TriangulationAndMeasure(Ω,v_cell_quad,is_nᵃ,is_a₋)
      Ωᶜ,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_nᶜ,is_c₋)

      dΩᶜ = Measure(Ωᶜ,2*order)
      nΓ = normal(φ₋,Ω)

      # Update buffer
      buffer[] = ( Ωᶜ=Ωᶜ,dΩᶜ=dΩᶜ,Ωˡ=Ωˡ,dΩˡ=dΩˡ,aggsˡ=aggsˡ,
                   dΓ=dΓ,nΓ=nΓ,cp₋=cp₋,φ₋=φ₋,t=t,Vbg=Vbg )
      return true

    end

  end

  # Reference FEs
  N = num_dims(bgmodel)
  reffeᵘ = ReferenceFE(lagrangian,VectorValue{N,Float64},order)
  reffeʷ = ReferenceFE(lagrangian,VectorValue{N,Float64},order-1)
  reffeˢ = ReferenceFE(lagrangian,VectorValue{N,Float64},order,space=:S)
  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1,space=:P)
  reffeᵉ = ReferenceFE(lagrangian,Float64,order-1)

  function update_all!(i::Int,t::Real,dt::Real,disp,val::Real)

    update_buffer!(i,t,dt,disp,val)

    # Triangulations
    Ωˡ = buffer[].Ωˡ
    Ωᶜ = buffer[].Ωᶜ

    # Aggregates
    aggsˡ = buffer[].aggsˡ

    # Measures and normal
    dΩˡ = buffer[].dΩˡ
    dΩᶜ = buffer[].dΩᶜ
    dΓ  = buffer[].dΓ
    nΓ  = buffer[].nΓ
    φ   = buffer[].φ₋

    # Ghost skeleton
    Λ  = SkeletonTriangulation(Ωᶜ)
    dΛ = Measure(Λ,2*order)
    nΛ = get_normal_vector(Λ)

    # Test FE spaces

    ## (u,p)-bulk
    Vstdᵘˡ = TestFESpace(Ωˡ,reffeᵘ,dirichlet_tags=["boundary"],
                                   dirichlet_masks=[(false,true)])
    Vserᵘˡ = TestFESpace(Ωˡ,reffeˢ,conformity=:L2)
    Vᵘˡ    = AgFEMSpace(Vstdᵘˡ,aggsˡ,Vserᵘˡ)
    Vstdᵖˡ = TestFESpace(Ωˡ,reffeᵖ)
    Vᵖˡ    = AgFEMSpace(Vstdᵖˡ,aggsˡ)

    # u-surface
    Vʷ = TestFESpace(Ωᶜ,reffeʷ,dirichlet_tags=["boundary"],
                               dirichlet_masks=[(false,true)])
    # e-surface
    Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
    # Lagrange multipliers
    Vˡ = ConstantFESpace(bgmodel)

    # Trial FE spaces
    Uᵘˡ = TrialFESpace(Vᵘˡ)
    Uᵖˡ = TrialFESpace(Vᵖˡ)
    Uʷ  = TrialFESpace(Vʷ)
    Uᵉ  = TrialFESpace(Vᵉ)
    Uˡ  = TrialFESpace(Vˡ)

    # Multifield FE spaces
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ])
    Yᵘ = MultiFieldFESpace([Vᵘˡ,Vᵖˡ,Vˡ])
    Xᵘ = MultiFieldFESpace([Uᵘˡ,Uᵖˡ,Uˡ])

    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ,dΛ,nΛ

  end

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0)
  m₀ = 2.0

  Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ,dΛ,nΛ = 
    update_all!(0,t₀,Δt,u₀,m₀)

  # *** WEAK FORM PARAMETERS ***
  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  # ** u-stabilisation **
  γʷ = γᶜ/h

  _eₕ = initial_density(Uᵉ,Xʳ,Yʳ,dΓ,dΩᶜ,nΓ)
  eₕ = interpolate_everywhere(_eₕ,Uᵉ)
  
  _uₕ(x) = VectorValue(0.0,0.0)
  _pₕ(x) = 0.0

  ulₕ = interpolate_everywhere(_uₕ,Xᵘ[1])
  plₕ = interpolate_everywhere(_pₕ,Xᵘ[2])
  υₕ  = interpolate_everywhere(_uₕ,Xᵛ[1])

  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  ps = PETScLinearSolver(mykspsetup)

  i = 0
  t = t₀
  
  tol = 1e-8

  while t < T + tol

    @info "Time step $i, time $t and time step $Δt"

    assemᵘ = SparseMatrixAssembler(Tm,Tv,Xᵘ,Yᵘ)
    assemᵛ = SparseMatrixAssembler(Tm,Tv,Xᵛ,Yᵛ)

    Aᵘ = nothing
    Aᵛ = nothing

    for i in 1:maxiter

      _υₕ = υₕ
      _ulₕ = ulₕ

      aᵛ,bᵛ = cortical_flow_problem_axisymmetric(
        ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,γʷ,Pe,μˡ,R,activity)
      Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
      υₕ,_ = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

      aᵘ,bᵘ = bulk_flow_problem_axisymmetric(
        υₕ,dΩˡ,dΓ,nΓ,μˡ,R,γᵅ,h)
      Aᵘ,Bᵘ = _assemble_problem(aᵘ,bᵘ,assemᵘ,Xᵘ,Yᵘ,Aᵘ)
      ulₕ,plₕ,_ = _solve_problem(Aᵘ,Bᵘ,Xᵘ,ps)

      chk1 = √( ∑( ∫( ((υₕ-_υₕ)⋅(υₕ-_υₕ))*y )dΓ ) ) / √( ∑( ∫( (_υₕ⋅_υₕ)*y )dΓ ) )
      chk2 = √( ∑( ∫( ((ulₕ-_ulₕ)⋅(ulₕ-_ulₕ))*y )dΓ ) ) / √( ∑( ∫( (_ulₕ⋅_ulₕ)*y )dΓ ) )
      chk = max(chk1,chk2)

      @info "chk1 = $chk1 and chk2 = $chk2 at i = $i"
      if chk < reltol
        break
      elseif i == maxiter
        error("No convergence")
      end

    end

    writesol && postprocess_all(φ,dΩˡ.quad.trian,dΩᶜ.quad.trian,
      eₕ,υₕ,ulₕ,plₕ,i=i,of=output_frequency,name=name)

    msₕ = get_maximum_magnitude_with_dirichlet(υₕ)

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ,dΛ,nΛ = 
      update_all!(i,t,Δt,υₕ,msₕ)

    # ** e-stabilisation **
    δ = 2.0 * msₕ * Δt
    γⁿ = γᶜ * ( msₕ + 1.0 / ( δ + h ) )
    γ¹ = γᶜ * msₕ * h
    @show msₕ, δ, γʷ, γⁿ, γ¹

    assemᵉ = SparseMatrixAssembler(Tm,Tv,Uᵉ,Vᵉ)
    aᵉ,bᵉ = transport_problem_axisymmetric(
      υₕ,eₕ,dΓ,nΓ,γⁿ,dΩᶜ,γ¹,dΛ,nΛ,τᵈkₒ,Δt)
    opᵉ = AffineFEOperator(aᵉ,bᵉ,Uᵉ,Vᵉ,assemᵉ)
    eₕ = solve(ps,opᵉ)

  end

end

struct RotatedEllipse
  c::VectorValue{2,Float64}
  a::Float64
  b::Float64
  θ::Float64
  φ::Function
  ∇φ::Function
  function RotatedEllipse(c::VectorValue{2,Float64},
                          a::Float64,
                          b::Float64,
                          θ::Float64)
    φ = ( x -> ( ( cos(θ) * ( x[1] - c[1] ) + sin(θ) * ( x[2] - c[2] ) ) / a )^2 + 
               ( ( cos(θ) * ( x[2] - c[2] ) - sin(θ) * ( x[1] - c[1] ) ) / b )^2 - 1.0 )
    ∇φ = x -> ∇(φ)(x)
    new(c,a,b,θ,φ,∇φ)
  end
end

function surface_bulk_viscous_flows_axisymmetric(
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            μˡ::Float64,
            R::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64,
            center::VectorValue{2,Float64},
            x_axis::Float64,
            y_axis::Float64;
            rotation::Float64 = 0.0,
            initial_density::Function = verification,
            activity::Function = unit_activity_axisymmetric,
            order::Int = 2,
            γᶜ::Float64 = 1.0,
            γᵅ::Float64 = 20.0,
            τᵈkₒ::Float64 = 10.0,
            writesol::Bool = true,
            reltol::Float64 = 1.0e-5,
            maxiter::Int = 10,
            output_frequency::Int = 1,
            redistance_frequency::Int = 1,
            name::String = "plt")

  # Background geometry for FE approximation
  h       = (domain[2]-domain[1])/n
  cells   = (n,div(n,2))
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω       = Triangulation(bgmodel)

  # Background geometry for Rigid Body (RB)
  domain  = (domain[1],domain[2],domain[1],domain[2])
  cells   = (n,n) # Symmetric domain
  rbmodel = CartesianDiscreteModel(domain,cells)
  Σ       = Triangulation(rbmodel)

  # Initial ellipse parameters
  a = x_axis; b = y_axis
  cˢ = center; θˢ = rotation
  cᵇ = center; θᵇ = rotation

  # Buffer of active model and integration objects
  degree = order < 3 ? 3 : 2*order
  buffer = Ref{Any}(( Ωᶜ  = nothing, dΩᶜ   = nothing,
                      Ωˡ  = nothing, dΩˡ   = nothing,
                      dΓ  = nothing, nΓ    = nothing,
                      φ₋  = nothing, aggsˡ = nothing,
                      cp₋ = nothing, t     = nothing,
                      dS  = nothing, dΣ    = nothing,
                      ψˢ  = nothing, ψᵇ    = nothing,
                      Vbg = nothing ))

  function update_buffer!(i,t,dt,v₋₂,mv₋₂)

    if buffer[].t == t
      return true
    else

      Ωᶜ    = buffer[].Ωᶜ
      dΩᶜ   = buffer[].dΩᶜ
      Ωˡ    = buffer[].Ωˡ
      dΩˡ   = buffer[].dΩˡ
      aggsˡ = buffer[].aggsˡ
      dΓ    = buffer[].dΓ
      nΓ    = buffer[].nΓ
      cp₋   = buffer[].cp₋
      φ₋    = buffer[].φ₋
      t     = buffer[].t
      Vbg   = buffer[].Vbg

      if buffer[].Ωᶜ === nothing
        Vbg  = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order))
        _φ₋  = interpolate_everywhere(ls.φ,Vbg)
      else
        cp₋₂ = buffer[].cp₋
        φ₋₂  = buffer[].φ₋
        __φ  = get_free_dof_values(φ₋₂.φ)
        Ωⱽ   = get_triangulation(Vbg)
        _ϕ₋  = compute_normal_displacement(cp₋₂,φ₋₂,v₋₂,dt,Ωⱽ)
        ϕ₋   = __φ - _ϕ₋
        _φ₋  = FEFunction(Vbg,ϕ₋)
      end

      # Current time level set
      φ₋ = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      ( i % redistance_frequency == 0 ) && begin
        _φ₋ = compute_distance_fe_function(bgmodel,Vbg,φ₋,order,cppdegree=3)
        φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      end

      cp₋ = compute_closest_point_projections(Vbg,φ₋,order,
              cppdegree=3,trim=true,limitstol=1.0e-2)

      # Current time surface and bulk measures
      squad = Quadrature(algoim,φ₋,degree,phase=CUT)
      s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)
      vquad = Quadrature(algoim,φ₋,degree,phase=IN)
      v_cell_quad,is_a₋ = CellQuadratureAndActiveMask(bgmodel,vquad)

      # Surface narrow-band triangulation
      δ₋ = 2.0 * mv₋₂ * dt
      _,is_nᶜ = narrow_band_triangulation(Ω,_φ₋,Vbg,is_c₋,δ₋)

      # Narrow band bulk active triangulations and measures
      # For bulk, only need to extend along positive LS vals
      _,is_nᵃ = active_triangulation(Ω,_φ₋,Vbg,is_a₋,δ₋)

      # Current aggregates
      aggsˡ = aggregate(Ω,is_a₋,is_c₋,IN,is_nᵃ)

      Ωˡ,dΩˡ = TriangulationAndMeasure(Ω,v_cell_quad,is_nᵃ,is_a₋)
      Ωᶜ,dΓ  = TriangulationAndMeasure(Ω,s_cell_quad,is_nᶜ,is_c₋)

      dΩᶜ = Measure(Ωᶜ,2*order)
      nΓ  = normal(φ₋,Ω)

      # Compute measures for rigid body
      _ψˢ   = RotatedEllipse(cˢ,a,b,θˢ)
      ψˢ    = AlgoimCallLevelSetFunction(_ψˢ.φ,_ψˢ.∇φ)
      squad = Quadrature(algoim,ψˢ,degree,phase=CUT)
      _,dS  = TriangulationAndMeasure(Σ,squad)

      _ψᵇ   = RotatedEllipse(cᵇ,a,b,θᵇ)
      ψᵇ    = AlgoimCallLevelSetFunction(_ψᵇ.φ,_ψᵇ.∇φ)
      vquad = Quadrature(algoim,ψᵇ,degree,phase=IN)
      _,dΣ  = TriangulationAndMeasure(Σ,vquad)

      # Update buffer
      buffer[] = ( Ωᶜ=Ωᶜ,dΩᶜ=dΩᶜ,Ωˡ=Ωˡ,dΩˡ=dΩˡ,aggsˡ=aggsˡ,dΓ=dΓ,
        nΓ=nΓ,cp₋=cp₋,φ₋=φ₋,t=t,dS=dS,dΣ=dΣ,ψˢ=ψˢ,ψᵇ=ψᵇ,Vbg=Vbg )
      return true

    end

  end

  # Reference FEs
  N = num_dims(bgmodel)
  reffeᵘ = ReferenceFE(lagrangian,VectorValue{N,Float64},order)
  reffeʷ = ReferenceFE(lagrangian,VectorValue{N,Float64},order-1)
  reffeˢ = ReferenceFE(lagrangian,VectorValue{N,Float64},order,space=:S)
  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1,space=:P)
  reffeᵉ = ReferenceFE(lagrangian,Float64,order-1)

  function update_all!(i::Int,t::Real,dt::Real,disp,val::Real)

    update_buffer!(i,t,dt,disp,val)

    # Triangulations and aggregates
    Ωˡ = buffer[].Ωˡ
    Ωᶜ = buffer[].Ωᶜ

    aggsˡ = buffer[].aggsˡ

    # Measures and normal
    dΩˡ = buffer[].dΩˡ
    dΩᶜ = buffer[].dΩᶜ
    dΓ  = buffer[].dΓ
    nΓ  = buffer[].nΓ
    φ   = buffer[].φ₋

    # Test FE spaces

    ## (u,p)-bulk
    Vstdᵘˡ = TestFESpace(Ωˡ,reffeᵘ,dirichlet_tags=["boundary"],
                                   dirichlet_masks=[(false,true)])
    Vserᵘˡ = TestFESpace(Ωˡ,reffeˢ,conformity=:L2)
    Vᵘˡ    = AgFEMSpace(Vstdᵘˡ,aggsˡ,Vserᵘˡ)
    Vstdᵖˡ = TestFESpace(Ωˡ,reffeᵖ)
    Vᵖˡ    = AgFEMSpace(Vstdᵖˡ,aggsˡ)

    # u-surface
    Vʷ = TestFESpace(Ωᶜ,reffeʷ,dirichlet_tags=[5],
                               dirichlet_masks=[(false,true)])
    # e-surface
    Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
    # Lagrange multipliers
    Vˡ = ConstantFESpace(bgmodel)

    # Trial FE spaces
    Uᵘˡ = TrialFESpace(Vᵘˡ)
    Uᵖˡ = TrialFESpace(Vᵖˡ)
    Uʷ  = TrialFESpace(Vʷ)
    Uᵉ  = TrialFESpace(Vᵉ)
    Uˡ  = TrialFESpace(Vˡ)

    # Multifield FE spaces
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ])
    Yᵘ = MultiFieldFESpace([Vᵘˡ,Vᵖˡ,Vˡ])
    Xᵘ = MultiFieldFESpace([Uᵘˡ,Uᵖˡ,Uˡ])

    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ

  end

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0)
  m₀ = 0.01

  Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ = update_all!(0,t₀,Δt,u₀,m₀)

  # *** WEAK FORM PARAMETERS ***
  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  # ** u-stabilisation **
  γʷ = γᶜ/h
  # ** e-stabilisation **
  γᵉ = γᶜ/h

  _eₕ = initial_density(Uᵉ,Xʳ,Yʳ,dΓ,dΩᶜ,nΓ)
  eₕ  = interpolate_everywhere(x->_eₕ(x,t₀),Uᵉ)
  
  _uₕ(x) = VectorValue(0.0,0.0)
  _pₕ(x) = 0.0

  ulₕ = interpolate_everywhere(_uₕ,Xᵘ[1])
  plₕ = interpolate_everywhere(_pₕ,Xᵘ[2])
  υₕ  = interpolate_everywhere(_uₕ,Xᵛ[1])

  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  ps = PETScLinearSolver(mykspsetup)

  i = 0
  t = t₀
  
  tol = 1e-8

  while t < T + tol

    @info "Time step $i, time $t and time step $Δt"

    assemᵘ = SparseMatrixAssembler(Tm,Tv,Xᵘ,Yᵘ)
    assemᵛ = SparseMatrixAssembler(Tm,Tv,Xᵛ,Yᵛ)

    Aᵘ = nothing
    Aᵛ = nothing

    for i in 1:maxiter

      _υₕ  = υₕ
      _ulₕ = ulₕ

      aᵛ,bᵛ = cortical_flow_problem_axisymmetric(
        ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,γʷ,Pe,μˡ,R,activity)
      Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
      υₕ,_  = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

      aᵘ,bᵘ = bulk_flow_problem_axisymmetric(
        υₕ,dΩˡ,dΓ,nΓ,μˡ,R,γᵅ,h)
      Aᵘ,Bᵘ     = _assemble_problem(aᵘ,bᵘ,assemᵘ,Xᵘ,Yᵘ,Aᵘ)
      ulₕ,plₕ,_ = _solve_problem(Aᵘ,Bᵘ,Xᵘ,ps)

      chk1 = √( ∑( ∫( ((υₕ-_υₕ)⋅(υₕ-_υₕ))*y )dΓ ) ) / √( ∑( ∫( (_υₕ⋅_υₕ)*y )dΓ ) )
      chk2 = √( ∑( ∫( ((ulₕ-_ulₕ)⋅(ulₕ-_ulₕ))*y )dΓ ) ) / √( ∑( ∫( (_ulₕ⋅_ulₕ)*y )dΓ ) )
      chk = max(chk1,chk2)

      @info "chk1 = $chk1 and chk2 = $chk2 at i = $i"
      if chk < reltol
        break
      elseif i == maxiter
        error("No convergence")
      end

    end

    msₕ = get_maximum_magnitude_with_dirichlet(υₕ)

    aυₕ = ( ∑( ∫( (ulₕ⋅VectorValue(1.0,0.0))*y )dΩˡ ) / ∑( ∫( y )dΩˡ ) ) * VectorValue(1.0,0.0)
     υₕ =  υₕ - aυₕ
    ulₕ = ulₕ - aυₕ
    @show aυₕ

    # Compute averages in the meridian
    dS = buffer[].dS
    dΣ = buffer[].dΣ

    Σˢ = dS.quad.trian
    Σᵇ = dΣ.quad.trian

    Vˢ = TestFESpace(Σˢ,reffeᵘ)
    Vᵇ = TestFESpace(Σᵇ,reffeᵘ)
    
    _uʳ(x) = x[2] > 0.0 ? ulₕ(x) : 
      VectorValue(ulₕ(Point(x[1],-x[2]))⋅VectorValue(1.0,0.0),
             -1.0*ulₕ(Point(x[1],-x[2]))⋅VectorValue(0.0,1.0))
     uˢ = interpolate_everywhere(x->_uʳ(x),Vˢ)
     uᵇ = interpolate_everywhere(x->_uʳ(x),Vᵇ)

    _θˢ = θˢ
    _θᵇ = θᵇ

    θˢ = θˢ + 0.5 * Δt * ( ∑( ∫( curl(uˢ) )dS ) / ∑( ∫( 1.0 )dS ) )
    θᵇ = θᵇ + 0.5 * Δt * ( ∑( ∫( curl(uᵇ) )dΣ ) / ∑( ∫( 1.0 )dΣ ) )

    νˢ = ∑( ∫( uˢ )dS ) / ∑( ∫( 1.0 )dS )
    νᵇ = ∑( ∫( uᵇ )dΣ ) / ∑( ∫( 1.0 )dΣ )

    cˢ = cˢ + Δt * νˢ
    cᵇ = cᵇ + Δt * νᵇ

    @show νˢ
    @show cˢ
    @show θˢ,θˢ-_θˢ
    
    @show νᵇ
    @show cᵇ
    @show θᵇ,θᵇ-_θᵇ

    writesol && postprocess_all(φ,dΩˡ.quad.trian,dΩᶜ.quad.trian,
      eₕ,υₕ,ulₕ,plₕ,i=i,of=output_frequency,name=name)

    ψˢ = buffer[].ψˢ
    writesol && if ( i % output_frequency == 0 )
      writevtk(Σˢ,name*"_rbs_$i",cellfields=["LS"=>ψˢ.φ,"curl"=>curl(uˢ)],nsubcells=4)
    end

    ψᵇ = buffer[].ψᵇ
    writesol && if ( i % output_frequency == 0 )
      writevtk(Σᵇ,name*"_rbb_$i",cellfields=["LS"=>ψᵇ.φ,"curl"=>curl(uᵇ)],nsubcells=4)
    end

    xΓ = dΓ.quad.cell_point.values
    xΓ = lazy_map(Reindex(xΓ),dΓ.quad.cell_point.ptrs)
    
    θΓ = lazy_map(Broadcasting(x->atan(x[2],x[1])),xΓ)
    θΓ = vcat(θΓ...)

    vΓ = lazy_map(υₕ,xΓ)
    eΓ = lazy_map(eₕ,xΓ)

    xΓ  = vcat(xΓ...)
    xxΓ = map(x->x.data[1],xΓ)
    xyΓ = map(x->x.data[2],xΓ)
  
    vΓ  = vcat(vΓ...)
    nvΓ = map(x->norm(x),vΓ)
    eΓ  = vcat(eΓ...)
     
    perm = sortperm(θΓ)
    θΓ  = θΓ[perm]
    xxΓ = xxΓ[perm]
    xyΓ = xyΓ[perm]
    nvΓ = nvΓ[perm]
    eΓ  = eΓ[perm]

    xxΓ = vcat(xxΓ,reverse(xxΓ))
    xyΓ = vcat(xyΓ,-1.0 .* reverse(xyΓ))
    nvΓ = vcat(nvΓ,-1.0 .* reverse(nvΓ))
    θΓ  = vcat(θΓ,2*π .- reverse(θΓ))
    eΓ  = vcat(eΓ,reverse(eΓ))
    
    cds = DataFrame("x"=>xxΓ,"y"=>xyΓ,"v"=>nvΓ,"theta"=>θΓ,"e"=>eΓ)
    CSV.write("examples/OocyteAndSpindle/results/results_$i.csv",cds)

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ = 
      update_all!(i,t,Δt,υₕ,msₕ)

    eₕ  = interpolate_everywhere(x->_eₕ(x,t),Uᵉ)

  end

end
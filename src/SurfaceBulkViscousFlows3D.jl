function surface_bulk_viscous_flows_3D(
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            μˡ::Float64,
            R::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64;
            initial_density::Function = verification,
            activity::Function = unit_activity,
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
  cells = (n,n,n)
  h = (domain[2]-domain[1])/n
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(bgmodel)

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
        __φ = get_free_dof_values(φ₋₂.φ)
        Ωⱽ  = get_triangulation(Vbg)
        _ϕ₋ = compute_normal_displacement(cp₋₂,φ₋₂,v₋₂,dt,Ωⱽ)
        ϕ₋  = __φ - _ϕ₋
        _φ₋ = FEFunction(Vbg,ϕ₋)
      end

      # Current time level set
      φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      ( i % redistance_frequency == 0 ) && begin
        _φ₋ = compute_distance_fe_function(bgmodel,Vbg,φ₋,order,cppdegree=3)
        φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      end

      cp₋ = compute_closest_point_projections(Vbg,φ₋,order,cppdegree=3)

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

    # Triangulations and aggregates
    Ωˡ = buffer[].Ωˡ
    Ωᶜ = buffer[].Ωᶜ

    aggsˡ    = buffer[].aggsˡ

    # Measures and normal
    dΩˡ = buffer[].dΩˡ
    dΩᶜ = buffer[].dΩᶜ
    dΓ  = buffer[].dΓ
    nΓ  = buffer[].nΓ
    φ   = buffer[].φ₋

    # Test FE spaces

    ## (u,p)-bulk
    Vstdᵘˡ = TestFESpace(Ωˡ,reffeᵘ)
    Vserᵘˡ = TestFESpace(Ωˡ,reffeˢ,conformity=:L2)
    Vᵘˡ = AgFEMSpace(Vstdᵘˡ,aggsˡ,Vserᵘˡ)
    Vstdᵖˡ = TestFESpace(Ωˡ,reffeᵖ)
    Vᵖˡ = AgFEMSpace(Vstdᵖˡ,aggsˡ)

    # u-surface
    Vʷ = TestFESpace(Ωᶜ,reffeʷ)
    # e-surface
    Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
    # Lagrange multipliers
    Vˡ = ConstantFESpace(bgmodel)

    # Trial FE spaces
    Uᵘˡ = TrialFESpace(Vᵘˡ)
    Uᵖˡ = TrialFESpace(Vᵖˡ)
    Uʷ = TrialFESpace(Vʷ)
    Uᵉ = TrialFESpace(Vᵉ)
    Uˡ = TrialFESpace(Vˡ)

    # Multifield FE spaces
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ,Vˡ,Vˡ,Vˡ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ,Uˡ,Uˡ,Uˡ,Uˡ,Uˡ])
    Yᵘ = MultiFieldFESpace([Vᵘˡ,Vᵖˡ,Vˡ])
    Xᵘ = MultiFieldFESpace([Uᵘˡ,Uᵖˡ,Uˡ])

    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ

  end

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0,0.0)
  m₀ = 2.0

  Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ = update_all!(0,t₀,Δt,u₀,m₀)

  # *** WEAK FORM PARAMETERS ***
  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  # ** u-stabilisation **
  γʷ = γᶜ/h
  # ** e-stabilisation **
  γᵉ = γᶜ/h

  _eₕ = initial_density(Uᵉ,Xʳ,Yʳ,dΓ,dΩᶜ,nΓ)
  eₕ = interpolate_everywhere(_eₕ,Uᵉ)
  
  _uₕ(x) = VectorValue(0.0,0.0,0.0)
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

      aᵛ,bᵛ = cortical_flow_problem_3D(
        ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,γʷ,Pe,μˡ,R,activity)
      Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
      υₕ,_ = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

      aᵘ,bᵘ = bulk_flow_problem_3D(υₕ,dΩˡ,dΓ,nΓ,μˡ,R,γᵅ,h)
      Aᵘ,Bᵘ = _assemble_problem(aᵘ,bᵘ,assemᵘ,Xᵘ,Yᵘ,Aᵘ)
      ulₕ,plₕ,_ = _solve_problem(Aᵘ,Bᵘ,Xᵘ,ps)

      chk1 = √( ∑( ∫( (υₕ-_υₕ)⋅(υₕ-_υₕ) )dΓ ) ) / √( ∑( ∫( _υₕ⋅_υₕ )dΓ ) )
      chk2 = √( ∑( ∫( (ulₕ-_ulₕ)⋅(ulₕ-_ulₕ) )dΓ ) ) / √( ∑( ∫( _ulₕ⋅_ulₕ )dΓ ) )
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

    msₕ = get_maximum_magnitude_3D(υₕ)

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ = 
      update_all!(i,t,Δt,υₕ,msₕ)

    assemᵉ = SparseMatrixAssembler(Tm,Tv,Uᵉ,Vᵉ)
    aᵉ,bᵉ = transport_problem_3D(υₕ,eₕ,dΓ,dΩᶜ,nΓ,Δt,γᵉ,τᵈkₒ)
    opᵉ = AffineFEOperator(aᵉ,bᵉ,Uᵉ,Vᵉ,assemᵉ)
    eₕ = solve(ps,opᵉ)

  end

end

struct RotatedEllipsoid
  p::VectorValue{3,Float64}
  a::Float64
  b::Float64
  c::Float64
  θ::Float64
  φ::Function
  ∇φ::Function
  function RotatedEllipsoid(p::VectorValue{3,Float64},
                            a::Float64,
                            b::Float64,
                            c::Float64,
                            θ::Float64)
    φ = ( x -> ( ( cos(θ) * ( x[1] - p[1] ) + sin(θ) * ( x[3] - p[3] ) ) / a )^2 +
               ( x[2] / b )^2 + 
               ( ( cos(θ) * ( x[3] - p[3] ) - sin(θ) * ( x[1] - p[1] ) ) / c )^2 - 1.0 )
    ∇φ = x -> ∇(φ)(x)
    new(p,a,b,c,θ,φ,∇φ)
  end
end

function surface_bulk_viscous_flows_3D(
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            μˡ::Float64,
            R::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64,
            center::VectorValue{3,Float64},
            x_axis::Float64,
            y_axis::Float64,
            z_axis::Float64;
            rotation::Float64 = 0.0,
            initial_density::Function = verification,
            activity::Function = unit_activity,
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
  h       = (domain[2]-domain[1])/n
  cells   = (n,n,n)
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω       = Triangulation(bgmodel)

  # Initial ellipse parameters
  a = x_axis; b = y_axis; c = z_axis
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
      φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      ( i % redistance_frequency == 0 ) && begin
        _φ₋ = compute_distance_fe_function(bgmodel,Vbg,φ₋,order,cppdegree=3)
        φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      end

      cp₋ = compute_closest_point_projections(Vbg,φ₋,order,cppdegree=3)

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
      _ψˢ   = RotatedEllipsoid(cˢ,a,b,c,θˢ)
      ψˢ    = AlgoimCallLevelSetFunction(_ψˢ.φ,_ψˢ.∇φ)
      squad = Quadrature(algoim,ψˢ,degree,phase=CUT)
      _,dS  = TriangulationAndMeasure(Ω,squad)

      _ψᵇ   = RotatedEllipsoid(cᵇ,a,b,c,θᵇ)
      ψᵇ    = AlgoimCallLevelSetFunction(_ψᵇ.φ,_ψᵇ.∇φ)
      vquad = Quadrature(algoim,ψᵇ,degree,phase=IN)
      _,dΣ  = TriangulationAndMeasure(Ω,vquad)

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
    Vstdᵘˡ = TestFESpace(Ωˡ,reffeᵘ)
    Vserᵘˡ = TestFESpace(Ωˡ,reffeˢ,conformity=:L2)
    Vᵘˡ    = AgFEMSpace(Vstdᵘˡ,aggsˡ,Vserᵘˡ)
    Vstdᵖˡ = TestFESpace(Ωˡ,reffeᵖ)
    Vᵖˡ    = AgFEMSpace(Vstdᵖˡ,aggsˡ)

    # u-surface
    Vʷ = TestFESpace(Ωᶜ,reffeʷ)
    # e-surface
    Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
    # Lagrange multipliers
    Vˡ = ConstantFESpace(bgmodel)

    # Trial FE spaces
    Uᵘˡ = TrialFESpace(Vᵘˡ)
    Uᵖˡ = TrialFESpace(Vᵖˡ)
    Uʷ = TrialFESpace(Vʷ)
    Uᵉ = TrialFESpace(Vᵉ)
    Uˡ = TrialFESpace(Vˡ)

    # Multifield FE spaces
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ,Vˡ,Vˡ,Vˡ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ,Uˡ,Uˡ,Uˡ,Uˡ,Uˡ])
    Yᵘ = MultiFieldFESpace([Vᵘˡ,Vᵖˡ,Vˡ])
    Xᵘ = MultiFieldFESpace([Uᵘˡ,Uᵖˡ,Uˡ])

    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ

  end

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0,0.0)
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
  
  _uₕ(x) = VectorValue(0.0,0.0,0.0)
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

      aᵛ,bᵛ = cortical_flow_problem_3D(
        ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,γʷ,Pe,μˡ,R,activity)
      Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
      υₕ,_ = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

      aᵘ,bᵘ = bulk_flow_problem_3D(υₕ,dΩˡ,dΓ,nΓ,μˡ,R,γᵅ,h)
      Aᵘ,Bᵘ = _assemble_problem(aᵘ,bᵘ,assemᵘ,Xᵘ,Yᵘ,Aᵘ)
      ulₕ,plₕ,_ = _solve_problem(Aᵘ,Bᵘ,Xᵘ,ps)

      chk1 = √( ∑( ∫( (υₕ-_υₕ)⋅(υₕ-_υₕ) )dΓ ) ) / √( ∑( ∫( _υₕ⋅_υₕ )dΓ ) )
      chk2 = √( ∑( ∫( (ulₕ-_ulₕ)⋅(ulₕ-_ulₕ) )dΓ ) ) / √( ∑( ∫( _ulₕ⋅_ulₕ )dΓ ) )
      chk = max(chk1,chk2)

      @info "chk1 = $chk1 and chk2 = $chk2 at i = $i"
      if chk < reltol
        break
      elseif i == maxiter
        error("No convergence")
      end

    end

    msₕ = get_maximum_magnitude_3D(υₕ)

    aυₕ = ( ∑( ∫( ulₕ⋅VectorValue(1.0,0.0,0.0) )dΩˡ ) / ∑( ∫( 1.0 )dΩˡ ) ) * VectorValue(1.0,0.0,0.0)
     υₕ =  υₕ - aυₕ
    ulₕ = ulₕ - aυₕ
    @show aυₕ

    # Compute averages
    dS = buffer[].dS
    dΣ = buffer[].dΣ

    _θˢ = θˢ
    _θᵇ = θᵇ

    θˢ = θˢ + 0.5 * Δt * ( ∑( ∫( curl(ulₕ) )dS ) / ∑( ∫( 1.0 )dS ) ) ⋅
          VectorValue(0.0,1.0,0.0)
    θᵇ = θᵇ + 0.5 * Δt * ( ∑( ∫( curl(ulₕ) )dΣ ) / ∑( ∫( 1.0 )dΣ ) ) ⋅
          VectorValue(0.0,1.0,0.0)

    νˢ = ∑( ∫( ulₕ )dS ) / ∑( ∫( 1.0 )dS )
    νᵇ = ∑( ∫( ulₕ )dΣ ) / ∑( ∫( 1.0 )dΣ )

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

    Σˢ = dS.quad.trian
    Σᵇ = dΣ.quad.trian

    ψˢ = buffer[].ψˢ
    writesol && if ( i % output_frequency == 0 )
      writevtk(Σˢ,name*"_rbs_$i",cellfields=["LS"=>ψˢ.φ,"curl"=>curl(ulₕ)],nsubcells=4)
    end

    ψᵇ = buffer[].ψᵇ
    writesol && if ( i % output_frequency == 0 )
      writevtk(Σᵇ,name*"_rbb_$i",cellfields=["LS"=>ψᵇ.φ,"curl"=>curl(ulₕ)],nsubcells=4)
    end

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xᵘ,Yᵘ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩˡ,dΩᶜ,dΓ,nΓ,φ = 
      update_all!(i,t,Δt,υₕ,msₕ)

    eₕ  = interpolate_everywhere(x->_eₕ(x,t),Uᵉ)

  end

end
function get_maximum_magnitude_with_dirichlet(υₕ::SingleFieldFEFunction)
  cv = get_cell_dof_values(υₕ)
  length(cv) == 0 && return 0.0
  max_array = lazy_map(cv) do ccv
    mv = 0.0
    reshape(ccv,:,2)
    for i in 1:size(ccv,1)
      _mv = ccv[i,:]⋅ccv[i,:]
      if _mv > mv
        mv = _mv
      end
    end
    √(mv)
  end
  maximum(max_array)
end

function get_maximum_magnitude_3D(υₕ::SingleFieldFEFunction)
  fv = reshape(get_free_dof_values(υₕ),3,:)
  mv = 0.0
  for j in 1:size(fv,2)
    _mv = fv[:,j]⋅fv[:,j]
    if _mv > mv
      mv = _mv
    end
  end
  √(mv)
end

function postprocess_all(φ,Ωˡ,Ωᶜ,eₕ,υₕ,ulₕ,plₕ;i=0,of=1,name="plt")
  if ( i % of == 0 )
    writevtk(Ωᶜ,name*"_sur_$i",cellfields=["LS"=>φ.φ,"uₕ"=>υₕ,"eₕ"=>eₕ],nsubcells=4)
    writevtk(Ωˡ,name*"_blk_$i",cellfields=["LS"=>φ.φ,"uₕ"=>ulₕ,"pₕ"=>plₕ],nsubcells=4)
  end
end
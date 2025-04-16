# For definitions, see D. Schöllhammer and T. P. Fries, "Kirchhoff-Love 
# shell theory based on tangential differential calculus", Computational 
# Mechanics (2019) 64:113-131

left_project(u,n) = u - n⊗(n⋅u) # u⋅P
right_project(u,n) = u - (u⋅n)⊗n # P⋅u

to_tangent_vector(u,n) = u - n*(n⋅u) # Only vectors
to_tangent_tensor(u,n) = left_project(right_project(u,n),n)

xₜ(v,n) = to_tangent_vector(v,n)
Aₜ(v,n) = to_tangent_tensor(v,n)

directional_gradient(u,n) = left_project(∇(u),n)
symmetric_directional_gradient(u,n) = symmetric_part(directional_gradient(u,n))

∇ᵈ(u,n) = directional_gradient(u,n)
εᵈ(u,n) = symmetric_directional_gradient(u,n)

covariant_gradient(u,n) = to_tangent_tensor(∇(u),n)
symmetric_covariant_gradient(u,n) = symmetric_part(covariant_gradient(u,n))

∇ᶜ(u,n) = covariant_gradient(u,n)
εᶜ(u,n) = symmetric_covariant_gradient(u,n)

directional_divergence(u,n) = ∇⋅u - (∇(u)⋅n)⋅n
divᵈ(u,n) = directional_divergence(u,n)
divᶜ(u,n) = directional_divergence(u,n)
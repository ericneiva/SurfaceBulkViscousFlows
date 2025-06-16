unit_activity_axisymmetric(x) = x[2]

unit_activity_3D(x) = 1.0

tmp(x) = abs(x[1]) < 3.0/90 ? max(x[2],3.0/90) : x[2]

function contractile_ring_axisymmetric(x)
  Am = 10.0; Aₒ = 1.0; ω = 0.1;
  ( Aₒ + ( Am - Aₒ ) * ( exp( -x[1]*x[1] / ( 2.0 * ω * ω ) ) ) ) * tmp(x)
end

function contractile_ring_3D(x)
  Am = 10.0; Aₒ = 1.0; ω = 0.1;
  ( Aₒ + ( Am - Aₒ ) * ( exp( -x[1]*x[1] / ( 2.0 * ω * ω ) ) ) )
end
unit_activity_axisymmetric(x) = x[2]

unit_activity_3D(x) = 1.0

function contractile_ring_axisymmetric(x)
  Am = 10.0; Aₒ = 1.0; ω = 0.1;
  ( Aₒ + ( Am - Aₒ ) * ( exp( -x[1]*x[1] / ( 2.0 * ω * ω ) ) ) ) * x[2]
end

function contractile_ring_3D(x)
  Am = 10.0; Aₒ = 1.0; ω = 0.1;
  ( Aₒ + ( Am - Aₒ ) * ( exp( -x[1]*x[1] / ( 2.0 * ω * ω ) ) ) )
end
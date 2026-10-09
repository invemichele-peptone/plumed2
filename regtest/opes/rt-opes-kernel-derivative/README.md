# Kernel derivatives consistent with the bias

Kernels are truncated at `r2 = rc2` (in units of sigma), and `c = exp(-rc2 / 2)`.
CUTOFF_SCHEME selects how:

- PLAIN (default, as in v2.9): the value is energy-shifted, `h * (exp(-r2 / 2) - c)`,
  but the derivative is that of the force-shifted kernel. Forces are not exactly
  minus the gradient of the bias.
- PR (plumed/plumed2#1463): energy-shifted value and its derivative,
  `-h * exp(-r2 / 2) * delta_i / sigma_i^2`, which jumps to zero at the cutoff.
- ALTERNATIVE: force-shifted value, `h * (exp(-r2 / 2) - c * (1 + (rc2 - r2) / 2))`,
  and its derivative, `-h * (exp(-r2 / 2) - c) * delta_i / sigma_i^2`; both go to
  zero continuously at the cutoff.

This test deposits several kernels and compares Cartesian forces against
central finite differences of the bias energy, with updates disabled during
the comparison, for the two consistent schemes. Both OPES_METAD and
OPES_METAD_EXPLORE are tested. A cutoff of 3 makes any inconsistency between
kernel value and derivative large enough to detect reliably.

Adapted from plumed/plumed2#1463.

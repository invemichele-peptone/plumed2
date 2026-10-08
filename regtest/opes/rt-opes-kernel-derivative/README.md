# Kernel derivatives consistent with the bias

Kernels are truncated at `r2 = rc2` (in units of sigma) with a force shift,
`h * (exp(-r2 / 2) - c * (1 + (rc2 - r2) / 2))` with `c = exp(-rc2 / 2)`, so that
both the kernel and its derivative, `-h * (exp(-r2 / 2) - c) * delta_i / sigma_i^2`,
go to zero continuously at the cutoff.

This test deposits several kernels and compares Cartesian forces against
central finite differences of the bias energy, with updates disabled during
the comparison. Both OPES_METAD and OPES_METAD_EXPLORE are tested. A cutoff
of 3 makes any inconsistency between kernel value and derivative large
enough to detect reliably.

Adapted from plumed/plumed2#1463, where the same test is used with an
energy-shifted kernel, `h * (exp(-r2 / 2) - c)`, whose derivative is instead
discontinuous at the cutoff.

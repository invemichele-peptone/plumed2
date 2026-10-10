# Exact neighbor list

With `NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL`, the neighbor list is rebuilt before any
kernel outside it can get within KERNEL_CUTOFF, and Zed is updated using all kernels.
The result must then be the same as without the neighbor list, up to round-off.

This test runs two PLUMED instances along the same random walk of two distance CVs, one
with the exact neighbor list and one without it, and compares bias and forces at every
step, for OPES_METAD and OPES_METAD_EXPLORE. Adaptive sigma gives kernels of different
widths, and most of them end up outside the list. The test fails with `NLIST_ZED=PLAIN`
(Zed updated using only the kernels in the list) and with `NLIST_REBUILD=PLAIN`
(kernels within the cutoff missed between rebuilds).

The second setting (`merges`) uses a COMPRESSION_THRESHOLD close to KERNEL_CUTOFF and a
small list (`NLIST_PARAMETERS=1.1,0.2`). After a merge, the merged kernel can then merge
again with a kernel outside the list, which happens twice with OPES_METAD. The test fails
if the recursive merges search only the kernels in the list.

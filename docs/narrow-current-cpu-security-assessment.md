# Current CPU CKKS security assessment

This is the historical low-noise assessment. The recommended retuned CPU
profile is described in [the new qualification](narrow-ll14-retuned-profile.md).

Parameter-Selection assessed the retained-cache N=16384 LL14 CPU profile at
TFHEpp revision `cb7e3fb087d0d01ff8ef10cd1b817c2c4df80b03`.
The current parameters fail the 128-bit classical screen under MATZOV.

| Exposed dense family | log2(q) | Integer noise RMS | BDGL16 bits | MATZOV bits |
|---|---:|---:|---:|---:|
| Native automorphism keys | 435 | 2.8283016878 | 128.96 | 126.76 |
| DD master | 433 | 2.8283016878 | 129.56 | 127.32 |
| Direct selectors | 440 | 102.0019083399 | 128.97 | 126.79 |

Each entry is the minimum of primal uSVP, primal BDD, dual and dual hybrid
estimates; primal BDD is weakest in all six cases. The estimator is pinned to
`3e48ef421ec256afddb3e7d2249a77eab6e9ba12`, uses unbounded generic samples and a
4096 maximum beta, and ran under SageMath 10.6. MATZOV is that estimator's
default cost model; the assessment script preserves BDGL16 as its historical
CLI default and selects MATZOV explicitly.

The implementation truncates continuous normal samples toward zero.
Configured sigma 3.2 and 102.4 therefore yield the integer RMS values above,
using the ideal-normal second moment. The estimator models variance-matched
discrete Gaussians; it does not model the exact error PMF. The table gives
model-dependent attack-cost estimates, not a demonstrated break or an exact
security level. Passing BDGL16 narrowly does not justify a general 128-bit claim.

The q56 sparse switch uses one signed coefficient per public block: 30 blocks
of width 529 and one of width 514. Its seven-attack uniform-sparse proxy gives
143.89 bits under BDGL16, but fixes 16 positive and 15 negative coefficients
and omits the public block structure. No full MATZOV uniform-sparse run was made.
The separate block/rotation sensitivities test weights 9 through 15 at fixed
zeta 14688 and the `saturated-5-prefix` layout under both cost models. They
reuse geometry-only posterior histograms (4000 samples per weight, seed 2119)
and recompute decoder probabilities at q56 and RMS 2.828301687834186.
The lowest tested sensitivity is 137.0519822924393 bits under MATZOV's
optimistic joint-bound iid95 scenario. An extended weight15 beta/d grid found
no cheaper tested point. This is a bounded heuristic search, not a security
lower bound. The iid95 limits describe the posterior model, not confidence
in security. Other weights, layouts, masks and zeta values remain unsearched;
actual decoder correlations, feasible dictionary speedups, exact error shape,
ring structure and KDM/correlated evaluation-key assumptions remain unqualified.

Retune the dense exposure with a security margin, then repeat numerical and
performance qualification. The present performance milestone remains a
measurement of these current parameters and does not establish paper security
or performance parity. No keys or inputs were generated, no private keys or
key payloads were read, and no CKKS parameters changed during this assessment.
Only saved public key headers were read to verify key metadata.

Evidence is committed in Parameter-Selection revision `cf8690607495291394b4ac7a743b1ecb74109f51`:

- [Full numerical assessment and reproducible command manifests](https://github.com/virtualsecureplatform/Parameter-Selection/blob/cf8690607495291394b4ac7a743b1ecb74109f51/python/estimates/narrow_ckks_current_cpu_security.json)
- Assessment JSON SHA-256: `9c88bc2a685b7fe3ce2fd64dfa6dbfd5df282412044322be5a447fbda0957c52`.
- All 16 launched estimator jobs terminated with exit status zero. Historical
  default BDGL16 results and the joint decoder's default-noise regression
  agree with saved controls. Joint decoder self-tests and invalid noise/model
  argument checks passed.
- Raw logs and full structured grids are retained locally at
  `/home/work-bleaker/kmatsuoka/ckks-diagnostics/current-cpu-security-20261017`;
  their hashes, commands, exit receipts, source hashes, generic reports and
  per-scenario structured minima are embedded in the committed JSON.

The older `narrow-dd-mask42-security-scope-review.json` describes a historical
profile; this assessment covers the current retained-cache CPU profile.

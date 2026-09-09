# Eight-bit integer scheme-switch multiplication

The integer path supports unsigned 8-bit operands and exact 16-bit products.
It uses the same keys and cryptographic parameters as [the four-bit path](Integer4Multiplication.md).
For input width b, encode x/M with M=2^(2b). Public post-multiplication scaling
is 256*M; extraction remains q=2^48 and scale=2^40. Integer cross terms vanish
on reduction to the output torus. No secret key is used during evaluation.

```cpp
TFHEEncryptInteger<8>(lhs, x, tfhe_key);
TFHEEncryptInteger<8>(rhs, y, tfhe_key);
TFHECKKSIntegerMultiply<8>(product, lhs, rhs, switching_key);
unsigned result = TFHEDecryptIntegerProduct<8>(product, tfhe_key);
```

The existing `TFHECKKSInteger4Key` and key generation function are shared across
widths. The original four-bit evaluation API remains available.

The Parameter-Selection security estimates are unchanged: minimum classical
security 131.54 bits across the evaluated attacks. With a 10-sigma tail cutoff,
the conditional integer-unit error bound is 0.439649 < 0.5, and the union failure
bound is 2^-46.09 for one key generation and 256 fresh-input products. This is
an analytical ideal-normal-sampler bound, not a probability inferred from tests.
See [the detailed bound](integer8-failure-bound.json).

Reproduce the full-dimension boundary/random sample:

```sh
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=25769803776 -- build/test/ckks/ckks_integer_multiplication --8bit
```

The sample comprises eight
boundary cases and 248 deterministically selected random operand pairs. Keys
and encryption randomness are fresh. This is not an exhaustive full-dimension
test of all 65,536 pairs. `--8bit-small` exhaustively checks the domain with
small, insecure test dimensions, for arithmetic regression only.

Validation: all 256 full-dimension sampled products passed. Key generation took
8.626 seconds and evaluation took 419.124 seconds (1.637 seconds/product). Peak
RSS was 14.050 GiB under the 24 GiB address-space cap, with no OOM. Maximum
observed error was 0.000194995 integer units. All 65,536 small-dimension pairs
and the four-bit regression also passed. See [the run record](integer8-validation.json).

The matching native-interface CLPX benchmark and parameter-screen limitations
are recorded in [the comparison](SchemeSwitchComparison.md).

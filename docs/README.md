# TFHEpp Documentation

TFHEpp is a pure C++20 implementation of the TFHE / CGGI fully homomorphic
encryption (FHE) scheme. It targets fast Boolean-gate evaluation and includes
modern TFHE extensions such as circuit bootstrapping, programmable bootstrapping
with multiple LUT outputs, subset/partial keys, and modified Chen-style packing
("annihilate" key switching).

This documentation lives in `docs/` and is meant to complement the root
`README.md`. If you are new to TFHEpp, start with **Getting Started**.

## Contents

- `GettingStarted.md` — prerequisites, build, and a 5‑minute example.
- `Build.md` — build system, CMake options, and integration patterns.
- `Parameters.md` — parameter sets and security/compile‑time selection.
- `API.md` — core types, keys, encryption/decryption, and gate APIs.
- `Advanced.md` — LUT bootstrapping, circuit bootstrapping, subset keys, packing,
  Double Decomposition, BFV++, and other advanced features.
- `CKKS.md` — experimental CKKS types, dense bootstrapping, product eval keys,
  and practical lvl6 validation commands.
- `Serialization.md` — exporting/importing keys and ciphertexts with cereal.

## Quick Links

- Main public header: `include/tfhe++.hpp`
- Parameter/type definitions: `include/params.hpp`, `include/params/`
- Boolean gates: `include/tfhe/gate.hpp`
- Gate bootstrapping / LUTs: `include/tfhe/gatebootstrapping.hpp`
- Circuit bootstrapping: `include/tfhe/circuitbootstrapping.hpp`
- CKKS: `include/ckks/ckks.hpp`, `test/ckks/ckks_bootstrap_validation.cpp`,
  `test/ckks/ckks_bootstrap_workflow.cpp`
- Tutorial programs: `tutorial/`
- Tests/examples: `test/`

// Reproduce the rejected37-bit internal profile: retained encrypted Tanh
// misses the23-bit gate and exits nonzero. No production default selects it.
// See docs/narrow-two-layer-frequency-experiment.json for the qualification.
#define TFHEPP_TWO_LAYER_BASIS_EXPERIMENT
#include "ckks_narrow_ll14_basis_ntt_fbs.cpp"

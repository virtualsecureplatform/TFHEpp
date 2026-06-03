#pragma once

#include "bfv/bfv++.hpp"
#include "clpx/bfv-clpx.hpp"
#include "tfhe/circuitbootstrapping.hpp"
#include "tfhe/cloudkey.hpp"
#include "tfhe/cmuxmem.hpp"
#include "tfhe/detwfa.hpp"
#include "tfhe/externs/cloudkey.hpp"
#include "tfhe/externs/detwfa.hpp"
#include "tfhe/externs/key.hpp"
#include "tfhe/externs/keyswitch.hpp"
#include "tfhe/externs/tlwe.hpp"
#include "tfhe/externs/trgsw.hpp"
#include "tfhe/externs/trlwe.hpp"
#include "tfhe/gate.hpp"
#include "tfhe/gatebootstrapping.hpp"
#include "tfhe/homdecomp.hpp"
#include "io-packet.hpp"
#include "tfhe/key.hpp"
#include "tfhe/keyswitch.hpp"
#include "tfhe/largelut.hpp"
#include "params.hpp"
#include "clpx/params/SS2CLPX.hpp"
#include "tfhe/tlwe.hpp"
#include "tfhe/trgsw.hpp"
#include "tfhe/trlwe.hpp"

#ifndef __clang__
// Because of some resons (may be clang bug?) this will gives linking error
// caused by mismatching name mangling.
#include "tfhe/externs/circuitbootstrapping.hpp"
#include "tfhe/externs/gate.hpp"
#include "tfhe/externs/gatebootstrapping.hpp"
#endif

// BFV SIMD slot operations
#include "bfv/bfv-slots.hpp"
#include "bfv/bfv-digitext.hpp"
#include "bfv/bfv-c2s.hpp"
#include "bfv/bfv-bootstrapping.hpp"
#include "bfv/bfv-gbfv.hpp"
#include "ckks/ckks.hpp"

// Application
#include "tfhe/aes.hpp"
#include "tfhe/ascon.hpp"

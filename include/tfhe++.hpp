#pragma once

#include "bfv++.hpp"
#include "bfv-clpx.hpp"
#include "circuitbootstrapping.hpp"
#include "cloudkey.hpp"
#include "cmuxmem.hpp"
#include "detwfa.hpp"
#include "externs/cloudkey.hpp"
#include "externs/detwfa.hpp"
#include "externs/key.hpp"
#include "externs/keyswitch.hpp"
#include "externs/tlwe.hpp"
#include "externs/trgsw.hpp"
#include "externs/trlwe.hpp"
#include "gate.hpp"
#include "gatebootstrapping.hpp"
#include "homdecomp.hpp"
#include "io-packet.hpp"
#include "key.hpp"
#include "keyswitch.hpp"
#include "largelut.hpp"
#include "params.hpp"
#include "params/SS2CLPX.hpp"
#include "tlwe.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"

#ifndef __clang__
// Because of some resons (may be clang bug?) this will gives linking error
// caused by mismatching name mangling.
#include "externs/circuitbootstrapping.hpp"
#include "externs/gate.hpp"
#include "externs/gatebootstrapping.hpp"
#endif

// BFV SIMD slot operations
#include "bfv-slots.hpp"
#include "bfv-digitext.hpp"
#include "bfv-c2s.hpp"
#include "bfv-bootstrapping.hpp"
#include "bfv-gbfv.hpp"
#include "ckks.hpp"

// Application
#include "aes.hpp"
#include "ascon.hpp"

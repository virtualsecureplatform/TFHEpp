#!/usr/bin/env python3
"""Parameter-Selection screen for the measured CLPX forward/multiplication path.

The installed reverse-switch model predates the measured implementation.
Consequently this script deliberately makes no end-to-end failure claim.
"""
import argparse
import dataclasses
import json
import math
from pathlib import Path
import sys

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--parameter-selection',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
a=parser.parse_args()
sys.path.insert(0,str(a.parameter_selection/'python'))
from noiseestimation import clpx as model
from noiseestimation.params import clpx as defaults
from scipy.special import log_ndtr

# Override the older Parameter-Selection preset with the measured C++ gadget.
class ForwardRing(defaults.SS2CLPXlvl2param):
    l=3
    lₐ=3
    ℬbit=13
    ℬₐbit=13
    ℬ=2**13
    ℬₐ=2**13
class ForwardPBS:
    domainP=defaults.lvlhalfparam
    targetP=ForwardRing
class Packing:
    t=defaults.SS2CLPXlvl22param.t
    basebit=defaults.SS2CLPXlvl22param.basebit
    domainP=ForwardRing
    targetP=ForwardRing

forward=model.estimate_tlwes_to_clpx(validbit=8,num_multi=4,shift=4,w=20,
    iksP=defaults.lvl1hparam,bkP=ForwardPBS,sskP=Packing)
mult=model.estimate_clpx_multiplication_depth(
    initial_coefficient_variance=forward.packed_variance,validbit=8,
    num_multi=4,shift=4,w=20,max_multiplications=1,P=ForwardRing,
    iksP=defaults.lvl1hparam,bkP=ForwardPBS,sskP=Packing,model='paper')
step=mult.steps[-1]
threshold=model.paper_clpx_fixed_noise_threshold(P=ForwardRing,w=20)
z=threshold/math.sqrt(step.digit_value_variance)
log2tail=1+float(log_ndtr(-z))/math.log(2)
report={
    'scope':'Forward switch and one CLPX multiplication noise screen only; not end-to-end correctness certification.',
    'parameters':{'ring_n':2048,'torus_bits':64,'forward_gadget_levels':3,'forward_gadget_bits':13,
                  'reverse_gadget_levels':4,'reverse_gadget_bits':10,'validbit':8,'productbit':16,
                  'num_multi':4,'shift':4,'w':20,'reverse_numdigit':5,'reverse_basebit':4},
    'forward_coefficient_log2_variance':math.log2(forward.packed_variance),
    'product_coefficient_log2_variance':math.log2(step.coefficient_variance),
    'paper_postconversion_log2_variance':math.log2(step.digit_value_variance),
    'paper_threshold_standard_deviations':z,
    'paper_scalar_log2_gaussian_tail':log2tail,
    'paper_16_outputs_256_products_log2_union_tail':min(0,log2tail+math.log2(16*256)),
    'end_to_end_failure_target_verified':False,
    'limitations':[
        'Gaussian variance screen and paper threshold; correlations and FFT rounding are not a formal bound.',
        'The local reverse estimator models the older temps0/temps1 circuit and two carry digits.',
        'The measured C++ reverse circuit uses previous_scaled-2*fid_stage2, rounded-digit PBS and one base-16 carry digit.',
        'The scalar paper tail omits semantic failures at intermediate PBS calls and cannot certify the whole path.'
    ]
}
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

"""Independently decode retained LL14 basis outputs; print aggregate errors only.

Set PYTHONPATH to Parameter-Selection/python/estimates. Arguments:
checkpoint.bin outputs.bin report.json validation.json
"""
import sys,json,struct,math
from pathlib import Path
from narrow_ckks_basis_checkpoint_validate import Reader
from narrow_ckks_prefix_arithmetic import PhaseDecoder,canonical_decode
checkpoint,outputs,reportfile,dest=map(Path,sys.argv[1:]);report=json.loads(reportfile.read_text())
with checkpoint.open('rb') as f:
 r=Reader(f);meta=r.words(9);n=meta[1]
 if meta!=(0x4e415252504b3031,n,85,41,379,49,34,100,448) or n!=16384:raise ValueError('profile mismatch')
 length=r.words(1)[0]
 if length!=4 or r.read(length)!=b'tanh':raise ValueError('function mismatch')
 secrets=[[v if v<1<<447 else v-(1<<448) for v in r.polynomial(n,7)] for _ in range(2)]
with outputs.open('rb') as f:
 r=Reader(f)
 meta=r.words(9)
 if meta[:6]!=(0x4c4c313446425331,n,448,85,44,48):raise ValueError('output metadata mismatch')
 output_q,output_scale,square_q=meta[6:]
 if (output_q,output_scale,square_q) not in ((93,44,49),(85,37,48)):raise ValueError('unsupported output profile')
 if (report['output_log_q'],report['output_log_delta'],report['post_product_log_q'])!=(output_q,output_scale,square_q):raise ValueError('output/report profile mismatch')
 length=r.words(1)[0]
 if length>8:raise ValueError('invalid function length')
 function=r.read(length).decode()
 if function!=report['function'] or function not in ('sigmoid','gelu','tanh'):raise ValueError('function report mismatch')
 values=struct.unpack(r.endian+str(n)+'d',r.read(n*8));slots=[complex(values[2*i],values[2*i+1]) for i in range(n//2)]
 if any(not math.isfinite(v.real) or not math.isfinite(v.imag) or abs(v.real)>4 or abs(v.imag)>4 for v in slots):raise ValueError('invalid domain')
 cts=[[r.polynomial(n,7) for _ in range(2)] for _ in range(4)]
 if f.read(1):raise ValueError('trailing output')
dense=PhaseDecoder(secrets[0]);sparse=PhaseDecoder(secrets[1])
fresh=canonical_decode(dense.phase(*cts[0],85),44)
coeff=[math.ldexp(v,-44) for v in sparse.phase(*cts[1],48)]
output=canonical_decode(dense.phase(*cts[2],output_q),output_scale);square=canonical_decode(dense.phase(*cts[3],square_q),output_scale)
def target(x):
 if function=='sigmoid':return 1/(1+math.exp(-x))
 if function=='gelu':return 0.5*x*(1+math.erf(x/math.sqrt(2)))
 return math.tanh(x)
expected=[target(v.real) for v in slots]
metrics={'fresh_input_max_error':max(abs(a-b) for a,b in zip(fresh,slots)),'s2c_max_error':max(max(abs(coeff[i]-slots[i].real),abs(coeff[i+n//2]-slots[i].imag)) for i in range(n//2)),'max_error':max(abs(a-b) for a,b in zip(output,expected)),'post_product_max_error':max(abs(a-b*b) for a,b in zip(square,expected))}
d={'function':function,'n':n,'output_log_q':output_q,'output_log_delta':output_scale,'post_product_log_q':square_q,'scope':'Independent integer phases for full dense input, structured S2C output, real functional refresh and remaining-level square. All coefficients/slots; no key values printed; no security claim.','phase_coefficients_checked':4*n,'decoded_output_slots_checked':3*n//2,'metrics':metrics,'accuracy_passed':metrics['max_error']<2**-23 and metrics['post_product_max_error']<report['product_tolerance']}
dest.write_text(json.dumps(d,indent=2)+'\n');print(function,'independent error',metrics['max_error'],'square',metrics['post_product_max_error'],flush=True);raise SystemExit(0 if d['accuracy_passed'] else 1)

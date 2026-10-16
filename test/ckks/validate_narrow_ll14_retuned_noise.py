"""Check actual saved-key noise with independent integer phase arithmetic.

Set PYTHONPATH to Parameter-Selection/python/estimates. Reads diagnostic secrets,
prints only aggregate residual metrics, and generates no keys or ciphertexts.
This checks sampled encryption residuals, not a Gaussian PMF or security proof.
"""
import argparse
import json
import math
from pathlib import Path
import struct
from narrow_ckks_basis_checkpoint_validate import Reader
from narrow_ckks_prefix_arithmetic import PhaseDecoder


def measure(decoder, ciphertext, expected, q, sigma):
    phase = decoder.phase(*ciphertext, q)
    errors = [(a-b+(1 << (q-1))) % (1 << q)-(1 << (q-1))
              for a, b in zip(phase, expected)]
    rms = math.sqrt(sum(e*e for e in errors)/len(errors))
    # Broad aggregate check on the saved diagnostic sample. This catches
    # wrong noise scales and metadata without asserting a distribution proof.
    return dict(coefficients_checked=len(errors), configured_sigma=sigma,
                rms=rms, mean=sum(errors)/len(errors),
                max_absolute=max(map(abs, errors)),
                passed=0.95*sigma < rms < 1.05*sigma)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    parser.add_argument('selectors', type=Path)
    parser.add_argument('dd_key', type=Path)
    parser.add_argument('native_keys', type=Path)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    n = 16384
    with args.checkpoint.open('rb') as f:
        r = Reader(f)
        if r.words(9) != (0x4e415252504b3031,n,85,41,379,49,34,100,448):
            raise ValueError('checkpoint profile mismatch')
        length = r.words(1)[0]
        if length != 4 or r.read(length) != b'tanh':
            raise ValueError('checkpoint function mismatch')
        secrets = [[v if v < 1 << 447 else v-(1 << 448)
                    for v in r.polynomial(n,7)] for _ in range(2)]
    decoder = PhaseDecoder(secrets[0])
    results = {}
    # Pick a head selector with known zero plaintext, using the sparse secret.
    position = next(i for i in range(529) if secrets[1][i])
    selector_index = 1 if position == 0 else 0
    with (args.selectors/'head-0.bin').open('rb') as f:
        r = Reader(f)
        if r.words(9) != (0x4e4152524f573033,n,440,57,42,529,448,0,0):
            raise ValueError('retuned selector metadata mismatch')
        f.seek(selector_index*2*n*7*8, 1)
        ciphertext = [r.polynomial(n,7) for _ in range(2)]
    results['selector'] = measure(decoder, ciphertext, [0]*n, 440, 32768)
    native_file = sorted(args.native_keys.glob('[0-9]*.bin'))[0]
    with native_file.open('rb') as f:
        r = Reader(f)
        meta = r.words(8)
        sigma_bits = struct.unpack('<Q',struct.pack('<d',4096.0))[0]
        if meta[:6] != (0x4e4155544f303031,n,448,435,32,14) or meta[7] != sigma_bits:
            raise ValueError('retuned native metadata mismatch')
        origin = r.read(64)
        f.seek(13*2*n*7*8,1)
        ciphertext = [r.polynomial(n,7) for _ in range(2)]
    transformed = [0]*n
    for i,value in enumerate(secrets[0]):
        j = i*meta[6] % (2*n)
        transformed[j % n] = value if j < n else -value
    results['native'] = measure(decoder, ciphertext, transformed, 435, 4096)
    with args.dd_key.open('rb') as f:
        r = Reader(f)
        if r.words(6) != (0x4e44444e30303032,n,433,20,20,sigma_bits):
            raise ValueError('retuned DD noise metadata mismatch')
        length = r.words(1)[0]
        if length != 64 or r.read(length) != origin:
            raise ValueError('DD/native origin mismatch')
        if r.words(7) != (0x434b4b5343444431,n,1,433,20,20,4):
            raise ValueError('compact DD metadata mismatch')
        f.seek(21*22*2*n*4,1)
        ciphertext = [[0]*n for _ in range(2)]
        for secondary in range(22):
            for component in range(2):
                digits = struct.unpack(r.endian+str(n)+'i',r.read(n*4))
                shift = (21-secondary)*20
                for i,value in enumerate(digits):
                    if not -(1 << 19) <= value < (1 << 19):
                        raise ValueError('invalid compact DD digit')
                    ciphertext[component][i] += value << shift
    square = [-v for v in decoder.phase(secrets[0],[0]*n,433)]
    results['dd'] = measure(decoder, ciphertext, square, 433, 4096)
    passed = all(x['passed'] for x in results.values())
    args.report.write_text(json.dumps(dict(scope=__doc__, results=results,
                                          passed=passed),indent=2)+'\n')
    for name,result in results.items():
        print(name,'sampled residual RMS',result['rms'],'passed',result['passed'],flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

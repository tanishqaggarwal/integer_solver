#!/usr/bin/env python3
"""verify_endtoend.py -- the machine-checkable END-TO-END verifier for the chain

    QUBO E=0  =(A)=>  comb/ladder constraints  =(B)=>  k*G=T  =(C)=>  EQUATIONS.txt OK

Input is a candidate ground state expressed as the decoded selector bits (synthetic
mode) or a full/partial x-assignment (live mode).  The verifier:

  (a) reconstructs the full x_0..x_N via the established propagation
      (synthetic: synth_circuit.reconstruct(bits); live: the assignment IS the
       reconstruction, missing wires default to 0 exactly as solve_lab/checker.py does),
  (b) runs the exact integer checker (solve_lab/checker.py machinery) on the equations,
  (c) prints  satisfied k/N  and  RESULT: OK  iff every equation is 0 in Z.

Modes
  --synth BITS [SEED]   build a fresh SAME-SHAPE planted instance, decode its bits,
                        reconstruct, and check.  Also runs a NEGATIVE test (flip one
                        bit -> must FAIL at the target comparison).  No live key needed.
  --live PARTIAL.json   run the checker on solve_lab/EQUATIONS.txt with PARTIAL.json.
  (no args)             synth demo at 16/20/24 bits + live validation on the record
                        partial best/new_instance_partial_39026.json.
"""
import sys, os, re, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
ANNEAL = os.path.abspath(os.path.join(HERE, '..', '..'))
LAB = os.path.abspath(os.path.join(ANNEAL, '..'))
sys.path.insert(0, LAB); sys.path.insert(0, ANNEAL)
import checker                      # solve_lab/checker.py
import synth.cert_reduction.synth_circuit as sc
from synth.gen import make

import tempfile
SCRATCH = os.environ.get('CERT_SCRATCH', tempfile.gettempdir())


def _check(eqs_path, v):
    """Return (ok, n, fails) using the exact-integer checker machinery."""
    codes, _ = checker.load_equations(eqs_path)
    fails = checker.evaluate_all(codes, v)
    return (len(fails) == 0, len(codes), fails)


def run_synth(bits, seed=3, verbose=True):
    inst = make(bits, seed=seed)
    art = sc.build(inst)
    meta, eqs = art['meta'], art['eqs']
    eqs_path = os.path.join(SCRATCH, f'synth_eq_{bits}_{seed}.txt')
    sc.write_equations(eqs, eqs_path)

    # (a) decode the ground state to selector bits, reconstruct full assignment
    truebits = [(inst.k >> i) & 1 for i in range(bits)]
    v = sc.reconstruct(meta, truebits, inst.curve)

    # (b)+(c) exact check
    ok, n, fails = _check(eqs_path, v)
    if verbose:
        print(f"[synth {bits}-bit seed={seed}] p={inst.curve.p}  planted k={inst.k}")
        print(f"  variables x_0..x_{meta['nvars']-1}   equations={n}")
        print(f"  reconstruct(256 bits) -> full assignment; checker: "
              f"satisfied {n-len(fails)}/{n}  "
              + ("RESULT: OK" if ok else f"RESULT: FAIL {fails[:5]}"))

    # NEGATIVE control: flip one selector bit -> accumulator misses T -> must FAIL
    flip = truebits[:]; flip[0] ^= 1
    vneg = sc.reconstruct(meta, flip, inst.curve)
    okn, n2, failsn = _check(eqs_path, vneg)
    if verbose:
        print(f"  negative control (1 bit flipped): satisfied {n2-len(failsn)}/{n2}  "
              + ("UNEXPECTED OK" if okn else "RESULT: FAIL (as required)"))
    return dict(bits=bits, seed=seed, n=n, ok=ok, neg_fail=(not okn), nvars=meta['nvars'])


def run_live(partial_path, eqs_path=None):
    eqs_path = eqs_path or os.path.join(LAB, '..', 'EQUATIONS.txt')
    with open(partial_path) as f:
        d = json.load(f)
    v = [0] * checker.NVARS
    for kk, val in d.items():
        idx = int(kk[2:]) if kk.startswith('x_') else int(kk)
        v[idx] = int(val)
    t0 = time.time()
    ok, n, fails = _check(eqs_path, v)
    print(f"[live] {os.path.relpath(partial_path, LAB)} vs EQUATIONS.txt")
    print(f"  satisfied {n-len(fails)}/{n}  ({len(fails)} failing)  "
          + ("RESULT: OK" if ok else "RESULT: FAIL") + f"   ({time.time()-t0:.1f}s)")
    if fails:
        print(f"  first failing line indices: {fails[:10]}")
    return dict(ok=ok, n=n, fails=len(fails))


def main():
    a = sys.argv[1:]
    if a and a[0] == '--synth':
        bits = int(a[1]); seed = int(a[2]) if len(a) > 2 else 3
        run_synth(bits, seed); return 0
    if a and a[0] == '--live':
        r = run_live(a[1])
        return 0 if r['ok'] else 1

    print("=" * 74)
    print("END-TO-END VERIFIER  (bits -> full assignment -> exact checker)")
    print("=" * 74)
    print("\n--- SYNTHETIC same-shape instances (planted key, full chain closes) ---")
    for bits in (16, 20, 24):
        run_synth(bits, seed=3)
    print("\n--- LIVE instance: validate the reconstruction PATH on the record partial ---")
    rp = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    if os.path.exists(rp):
        run_live(rp)
        print("  NOTE: 7 equations remain open — the 256-bit ECDLP core is unsolved.")
        print("        Full OK on the LIVE instance requires the unknown selector bits;")
        print("        it is DEMONSTRATED above on same-shape synthetic instances.")
    return 0


if __name__ == '__main__':
    sys.exit(main())

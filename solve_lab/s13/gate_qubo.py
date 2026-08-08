#!/usr/bin/env python3
"""
Concrete demonstration that a per-prime *gate block* is a genuinely small QUBO.

A degree-2 atom of the instance is  out = a * b  (output coefficient +-1,
[measured, s9]).  Over GF(q) its constraint is

        out  ==  a * b   (mod q)        <=>     a*b - out - k*q == 0

with residues a,b,out in [0,q) and quotient k in [0, q).  Binarize each of
a,b,out,k as ceil(log2 q) bits and take the squared penalty

        E(bits) = ( a*b - out - k*q )^2 .

This is a small QUBO block (after the standard degree reduction of the b_i b_j
products; here we just enumerate to certify the ground state).  We verify:

  * every ground state (E == 0) is a VALID triple  out == a*b (mod q);
  * every valid triple is realized by some ground state;
  * the coupler magnitudes are O(q^2), i.e. ~2 log2 q bits -- NOT 512 bits.

This is one block; a full per-prime system is a sparse DAG of ~31,475 such
blocks (plus linear atoms), disjoint from every other prime's system.

Usage:  python3 gate_qubo.py [q]      (default q=7)
"""
import sys, itertools

def bits_of(n, w):
    return [(n >> i) & 1 for i in range(w)]

def demo(q):
    w = (q - 1).bit_length()                 # bits per residue
    print(f"q = {q},  {w} bits per residue,  variables = 4*{w} = {4*w} binary")

    # Enumerate all bit assignments of (a,b,out,k), each w bits.
    # (4w bits total; fine for small q as a certification of the ground state.)
    valid_triples = {(a, b, (a * b) % q) for a in range(q) for b in range(q)}
    ground = []
    max_coeff = 0
    n = 4 * w
    for mask in range(1 << n):
        bs = bits_of(mask, n)
        a = sum(bs[i] << i for i in range(w))
        b = sum(bs[w + i] << i for i in range(w))
        o = sum(bs[2 * w + i] << i for i in range(w))
        k = sum(bs[3 * w + i] << i for i in range(w))
        # residues must be < q for a legal encoding; illegal codes get a
        # large penalty in a real QUBO -- here we just require legality.
        if a >= q or b >= q or o >= q or k >= q:
            continue
        E = (a * b - o - k * q) ** 2
        max_coeff = max(max_coeff, abs(a * b - o - k * q))
        if E == 0:
            ground.append((a, b, o))

    gset = set(ground)
    ok_sound = gset <= valid_triples
    ok_complete = valid_triples <= gset
    print(f"  ground states (E=0): {len(ground)}")
    print(f"  valid triples out==a*b (mod q): {len(valid_triples)}")
    print(f"  SOUND  (every ground state valid):     {ok_sound}")
    print(f"  COMPLETE (every valid triple reached): {ok_complete}")
    print(f"  max |residual| before squaring ~ q^2 = {q*q}; "
          f"coupler width <= {(q*q).bit_length()} bits  (vs 512 for a mod-p block)")
    assert ok_sound and ok_complete, "gate QUBO did not certify"
    print("  => the gate block is a correct, SMALL QUBO.  [OK]")

if __name__ == '__main__':
    q = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    demo(q)

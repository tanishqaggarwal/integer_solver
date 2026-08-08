#!/usr/bin/env python3
"""bench.py -- does the encoding anneal? largest modmul / largest mu any solver cracks."""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath('.'))
from synth.gen import make
import model as M, solvers as S
import numpy as np

def hit(ising, e): return e == 0

def run_modmul(sizes, reps=5, effort='mid'):
    budg = dict(low=dict(sa=2000, pt=(1000,8), tabu=8000, sb=1500),
                mid=dict(sa=20000, pt=(6000,12), tabu=40000, sb=8000),
                hi =dict(sa=120000, pt=(30000,16), tabu=200000, sb=40000))[effort]
    print(f"\n== MODMUL (free search), effort={effort}, {reps} reps ==")
    print(f"{'s':>3} {'n':>6} {'sa':>10} {'pt':>10} {'tabu':>10} {'sb':>10}")
    out = {}
    for s in sizes:
        mm = M.build_modmul(s, mode='wallace')
        row = {}
        for name in ('sa','pt','tabu','sb'):
            best = 10**9; hits = 0
            for r in range(reps):
                if name=='sa': e,_=S.sa(mm['ising'], sweeps=budg['sa'], seed=r)
                elif name=='pt': e,_=S.pt(mm['ising'], sweeps=budg['pt'][0], R=budg['pt'][1], seed=r)
                elif name=='tabu': e,_=S.tabu(mm['ising'], iters=budg['tabu'], seed=r)
                else: e,_=S.sb(mm['ising'], steps=budg['sb'], seed=r)
                best=min(best,e); hits+= (e==0)
            row[name]=(hits,best)
        out[s]=row
        print(f"{s:3d} {mm['Q'].n:6d} " + " ".join(f"{row[k][0]}/{reps}:{row[k][1]:<4.0f}" for k in ('sa','pt','tabu','sb')))
    return out

def run_comb_clamped(bits, mus, w=2, reps=4, effort='mid'):
    """interval-split sub-instance: answer digits CLAMPED to planted (fill ancillas).
    This is the easiest possible task -- a unique forced completion, no search over k."""
    budg = dict(low=8000, mid=40000, hi=200000)[effort]
    print(f"\n== COMB, answer-CLAMPED (ancilla fill only), tabu iters={budg}, {reps} reps ==")
    print(f"{'bits':>4} {'mu':>3} {'n':>7} {'clamped E=0':>12} {'best E':>7}")
    inst = make(bits, seed=3)
    out = {}
    for mu in mus:
        md = M.build_comb(inst, mu, w=w, mode='wallace')
        import random
        hits=0; best=10**9
        for r in range(reps):
            x0, clamp = M.onehot_start(md, rng=random.Random(r))
            e,_ = S.tabu(md['ising'], iters=budg, seed=r, x0=x0, clamp=clamp)
            best=min(best,e); hits+=(e==0)
        out[mu]=(hits,best)
        print(f"{bits:4d} {mu:3d} {md['Q'].n:7d} {hits}/{reps:<11d} {best:7.0f}")
    return out

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv)>1 else 'modmul'
    eff = sys.argv[2] if len(sys.argv)>2 else 'mid'
    if which=='modmul':
        r=run_modmul([6,8,10,12,14,16], effort=eff)
        json.dump({str(k):v for k,v in r.items()}, open('synth/solver/modmul_bench.json','w'), indent=1)
    elif which=='comb':
        r=run_comb_clamped(16, [2,4,6,8,10,12], effort=eff)
        json.dump({str(k):v for k,v in r.items()}, open('synth/solver/comb_bench.json','w'), indent=1)
    print("done")

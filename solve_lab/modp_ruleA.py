#!/usr/bin/env python3
"""Mod-P propagation with Rule A (boolean forcing): for an atom with a boolean
unknown b, if setting b=0 makes it a nonzero constant then force b=1 (and vice
versa). Seeded with one bit=1, the residue-load may cascade to force a whole
consistent bit set. Test every seed bit; report violations."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars
from repair import boolean_vars
from modp import P, inv, substitute_modp, solve_single_modp, NVARS

def subst_one_modp(poly, var, value):
    out = defaultdict(int)
    for m, c in poly.items():
        coef = c % P; newm = []
        for v in m:
            if v == var: coef = (coef * value) % P
            else: newm.append(v)
        out[tuple(newm)] = (out[tuple(newm)] + coef) % P
    return {m: c for m, c in out.items() if c % P != 0}

class RuleAEngine:
    def __init__(self, atoms, bset):
        self.atoms = atoms; self.bset = bset
        self.var_atoms = defaultdict(list)
        for ai, poly in enumerate(atoms):
            for v in atom_vars(poly): self.var_atoms[v].append(ai)
        self.val = [None]*NVARS; self.contra = 0
        self.wl = deque(range(len(atoms))); self.inwl = [True]*len(atoms)
    def assign(self, v, x):
        x %= P
        if self.val[v] is not None:
            if self.val[v] != x: self.contra += 1
            return
        self.val[v] = x
        for ai in self.var_atoms[v]:
            if not self.inwl[ai]: self.inwl[ai]=True; self.wl.append(ai)
    def propagate(self, ruleA=True):
        while self.wl:
            ai = self.wl.popleft(); self.inwl[ai]=False
            poly = substitute_modp(self.atoms[ai], self.val)
            uv = atom_vars(poly)
            if len(uv)==0:
                if poly.get((),0)%P!=0: self.contra+=1
                continue
            if len(uv)==1:
                k,d = solve_single_modp(poly)
                if k=='val': self.assign(*d)
                elif k=='dom' and len(d[1])==1: self.assign(d[0], next(iter(d[1])))
                continue
            if ruleA and len(uv)<=4:
                for b in [x for x in uv if x in self.bset]:
                    p0 = subst_one_modp(poly, b, 0); p1 = subst_one_modp(poly, b, 1)
                    bad0 = (not atom_vars(p0)) and p0.get((),0)%P!=0
                    bad1 = (not atom_vars(p1)) and p1.get((),0)%P!=0
                    if bad0 and not bad1: self.assign(b,1); break
                    if bad1 and not bad0: self.assign(b,0); break
    def viol(self, matoms):
        val=[x if x is not None else 0 for x in self.val]; vi=0
        for poly in matoms:
            s=0
            for m,c in poly.items():
                t=c%P
                for x in m: t=(t*val[x])%P
                s=(s+t)%P
            if s%P!=0: vi+=1
        return vi

def main():
    atoms = load_atoms(); bset = boolean_vars(atoms)
    mainv = set(json.load(open('main_comp.json'))['main_vars'])
    matoms = [poly for poly in atoms if set().union(*[set(m) for m in poly])&mainv]
    control = json.load(open('control_bits.json'))
    t0=time.time()
    # baseline with Rule A, no seed
    e=RuleAEngine(atoms,bset); e.propagate()
    for v in [b for b in bset if e.val[b] is None]:
        if e.val[v] is None: e.assign(v,0); e.propagate()
    for v in range(NVARS):
        if e.val[v] is None: e.assign(v,0); e.propagate()
    print(f"Rule-A baseline (no seed): {e.viol(matoms)} viol, {e.contra} contra ({time.time()-t0:.0f}s)", flush=True)
    # base state after pins for seeding
    base=RuleAEngine(atoms,bset); base.propagate()
    bval=list(base.val)
    best=(999,None)
    for k,seed in enumerate(control):
        eng=RuleAEngine(atoms,bset); eng.val=list(bval); eng.wl=deque(); eng.inwl=[False]*len(atoms)
        eng.assign(seed,1); eng.propagate(ruleA=True)
        nforced=sum(1 for b in control if eng.val[b] is not None)
        for v in [b for b in bset if eng.val[b] is None]:
            if eng.val[v] is None: eng.assign(v,0); eng.propagate(ruleA=True)
        for v in range(NVARS):
            if eng.val[v] is None: eng.assign(v,0)
        eng.propagate(ruleA=False)
        vi=eng.viol(matoms)
        if vi<best[0]: best=(vi,seed); print(f"  seed x_{seed}: viol={vi}, forced {nforced} bits, contra {eng.contra} ({time.time()-t0:.0f}s)", flush=True)
        if vi==0:
            print("SOLVED via seeded Rule-A!");
            val=[x if x is not None else 0 for x in eng.val]
            json.dump({f"x_{i}":int(val[i]) for i in range(NVARS)}, open('cand_ruleA.json','w')); return
    print(f"best seeded Rule-A: {best}", flush=True)

if __name__=='__main__':
    main()

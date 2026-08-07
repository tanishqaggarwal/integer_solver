"""Whole-instance pricing of the TWO congruences.

  congruence 1:  A0 + 7376877*A6 == Ca  (mod 7376877*P),  Ca = x_7068 - x_2099
  congruence 2:  A1              == Cb  (mod P),          Cb = x_4432 - x_19964

Removing either one lowers c by 1 and hence the failing count by 1 (12 - 7 + c).
So: is there ANY variable in the instance whose perturbation moves Ca or Cb mod P at a
cost of 0 extra failing equations?  Scan all 38,748 variables.
"""
import os, sys, json, time, random, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
P=2**256-2**32-977
M=7376877*P
BROKEN_GATES={22229,22230,35758,35761,35762}
BASE={22229,22230,35758,35759,35760,35761,35762}
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
base_av=L.all_atom_values(v0)
base_fail=set(L.failing_eqs(base_av))
Ca0=(v0[7068]-v0[2099])%M
Cb0=(v0[4432]-v0[19964])%P
print('Ca0 mod M =',Ca0%P,' Cb0 =',Cb0)
random.seed(23)
DS=[1,random.randrange(1,2**64)]
lo,hi=int(sys.argv[1]),int(sys.argv[2])
out=[]
t0=time.time()
for t in range(lo,min(hi,L.NVARS)):
    for D in DS:
        w=list(v0)
        blk=set(BROKEN_GATES)
        d=L.definer.get(t)
        if d is not None: blk.add(d)
        try:
            ch,_=L.ripple(w,{t:v0[t]+D},maxsteps=60000,block=blk)
        except Exception: continue
        dCa=((w[7068]-w[2099])%M)!=Ca0
        dCb=((w[4432]-w[19964])%P)!=Cb0
        if not (dCa or dCb): continue
        cand=set()
        for u in ch: cand.update(L.var_atoms[u])
        nz=set(BASE); diff=set()
        for a in cand:
            nv=L.evalpoly(L.polys[a],w)
            if nv!=base_av[a]: diff.add(a)
            if nv: nz.add(a)
            else: nz.discard(a)
        eqs=set(base_fail)
        for a in diff: eqs.update(L.atom2eq.get(a,()))
        avc={}; fail=0
        for e in eqs:
            m,sq,co=L.eq_atoms[e]
            s=0
            for a,c in co.items():
                if a not in avc: avc[a]=L.evalpoly(L.polys[a],w)
                s+=c*avc[a]
            if s: fail+=1
        E=set()
        for a in nz: E.update(L.atom2eq.get(a,()))
        out.append((fail,len(E),len(nz),t,str(D)[:4],int(dCa),int(dCb),sorted(nz-BASE)[:8]))
    if (t-lo)%2000==0: print(' ',t,len(out),f'{time.time()-t0:.0f}s',flush=True)
out.sort()
print(f'\nperturbations that MOVE Ca or Cb: {len(out)}')
print('cheapest 30 (fail, |E|, |S|, var, delta, dCa, dCb, extra atoms):')
for r in out[:30]: print('  ',r)
json.dump(out[:500],open(os.path.join(HERE,f'pa_price_{lo}.json'),'w'))
print(f'{time.time()-t0:.0f}s done')

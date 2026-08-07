"""Exhaustive single-move argument for c >= 2.

A congruence can only relax if one of these three mod-p residues moves:
  R_A  = x_7068 - x_2099            (congruence 1's C0)
  R_B  = x_14853 - x_1308           (a29539's pin -> licenses x_7068 mod p)
  R_C  = x_24548 - x_25442          (a7930's pin  -> licenses x_28730 mod p)
Enumerate ALL free inputs in the ancestor cones of x_2099, x_1308, x_14853,
x_24548, x_25442 (in frame 2), and price every one of them exactly.
"""
import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]; SS=set(SEVEN)
E12=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds=2):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8)
def cone(u):
    seen=set(); st=[u]
    while st:
        t=st.pop()
        if t in seen: continue
        seen.add(t)
        a=definer.get(t)
        if a is None: continue
        for x in L.avars[a]:
            if x!=t and x not in seen: st.append(x)
    return seen
tg={'x_2099':2099,'x_1308':1308,'x_14853':14853,'x_24548':24548,'x_25442':25442}
cand=set()
for k,t in tg.items():
    c=cone(t); fr=set(u for u in c if u not in definer)
    print(f'{k}: cone {len(c)} vars, {len(fr)} free inputs -> {sorted(fr)}')
    cand|=fr
cand-= set(DETACH)          # detached ones are handled separately
print(f'\nTOTAL distinct free inputs that can reach any pin residue: {len(cand)}')
print(sorted(cand))
print('\nexact price + does it move a pin residue?  (frame 2, delta = 1000003)')
R=lambda v:( (v[7068]-v[2099])%P, (v[14853]-v[1308])%P, (v[24548]-v[25442])%P )
R0=R(w)
rows=[]
for u in sorted(cand):
    w2=list(w); w2[u]+=1000003; fwd2(w2,2)
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    mv=[n for n,(a,b) in zip(('RA','RB','RC'),zip(R(w2),R0)) if a!=b]
    rows.append((len(f-E12),u,mv,L.NEQ-len(f)))
    print(f'  x_{u:<6} out12={len(f-E12):<4} score={L.NEQ-len(f):<6} moves={mv}', flush=True)
rows.sort()
print('\nsorted by outside-12 cost:')
for c,u,mv,s in rows: print(f'  x_{u:<6} cost={c:<4} moves={mv}')
movers=[r for r in rows if r[2]]
print(f'\nfree inputs that MOVE a pin residue: {len(movers)}; cheapest = '
      f'{min(movers)[0] if movers else None} equations')
zero=[r for r in rows if r[0]==0]
print(f'zero-collateral among them: {[(r[1],r[2]) for r in zero]}')

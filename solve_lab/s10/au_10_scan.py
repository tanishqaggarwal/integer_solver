import os, sys, time, json, collections, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
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
FREE=[u for u in range(L.NVARS) if u not in definer]
b2099,b7068,b28730 = w[2099]%P, w[7068]%P, w[28730]%P
avb=L.all_atom_values(w)
res={}
t0=time.time()
DELT = 1000003
for i,u in enumerate(FREE):
    w2=list(w); w2[u]+=DELT; fwd2(w2,2)
    av=L.all_atom_values(w2); f=set(L.failing_eqs(av))
    out=sorted(f-E)
    moved=[a for a in SEVEN if av[a]!=avb[a]]
    rec={'out':len(out),'outeqs':out[:40],
         'd2099':int((w2[2099]-w[2099])%P!=0),
         'd7068':int((w2[7068]-w[7068])%P!=0),
         'd28730':int((w2[28730]-w[28730])%P!=0),
         'movedSEVEN':moved,'score':L.NEQ-len(f)}
    res[u]=rec
    if i%500==0:
        print(f'{i}/{len(FREE)} {time.time()-t0:.0f}s', flush=True)
json.dump(res, open(os.path.join(HERE,'au_scan.json'),'w'))
zero=[u for u in FREE if res[u]['out']==0]
print('ZERO-collateral free params:', len(zero))
for u in zero:
    r=res[u]
    print(f'  x_{u:<6} score={r["score"]} d2099={r["d2099"]} d7068={r["d7068"]} d28730={r["d28730"]} movedSEVEN={r["movedSEVEN"]}')
print('done %.0fs'%(time.time()-t0))

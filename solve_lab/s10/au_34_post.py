"""Post-process au_scan.json: for every zero-collateral free input, check whether it
moves ANY variable by an amount that is NOT a multiple of p.  If none does, the two
mod-p congruences cannot be relaxed by any single-parameter move."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]; SS=set(SEVEN)
E=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8)
res=json.load(open(os.path.join(HERE,'au_scan.json')))
zero=sorted(int(u) for u,r in res.items() if r['out']==0)
print('zero-collateral free params:', len(zero))
print(zero)
DIRECT={7068,28730,29854,31864,642}   # the detached ones -- expected
HANDLES={1329,9118,10903,8731,9413,17325}
print('\nof these: detached', sorted(set(zero)&DIRECT), ' known handles', sorted(set(zero)&HANDLES))
print('others:', sorted(set(zero)-DIRECT-HANDLES))
print('\nchecking whether each zero-collateral knob moves any variable by a NON-multiple of p')
bad=[]
for u in zero:
    w2=list(w); w2[u]+=1000003; fwd2(w2,3)
    nonp=[t for t in range(L.NVARS) if w2[t]!=w[t] and (w2[t]-w[t])%P!=0 and t!=u]
    r=res[str(u)]
    tag=''
    if u in DIRECT: tag=' [detached]'
    elif u in HANDLES: tag=' [handle]'
    print(f'  x_{u:<6}{tag} moved_nonmult_p={len(nonp)} {nonp[:12]}  movedSEVEN={r["movedSEVEN"]}', flush=True)
    if nonp and u not in DIRECT and u not in HANDLES: bad.append((u,nonp))
print('\nSUSPECT knobs (zero collateral, move something non-mult-p, not a known handle):', bad[:20])
# distribution of collateral cost
cnt=collections.Counter(r['out'] for r in res.values())
print('\ndistribution of #failing-equations-outside-the-12 over all free params:')
for k in sorted(cnt)[:15]: print(f'   out={k}: {cnt[k]}')
cheap=sorted(((r['out'],int(u)) for u,r in res.items() if 0<r['out']<=6))
print('\ncheapest non-zero-collateral knobs:', cheap[:25])
c0movers=sorted(int(u) for u,r in res.items() if r['d2099'])
print('\nfree params that move x_2099 mod p (i.e. move C0):')
for u in c0movers: print(f'   x_{u}: out={res[str(u)]["out"]}')

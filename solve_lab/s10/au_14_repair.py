import os, sys, collections
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
w=list(v0); fwd2(w,8); avb=L.all_atom_values(w)
def st(u):
    d=L.definer.get(u); return f'x_{u} free={d is None} def=a{d} val={str(w[u])[:36]} bits={w[u].bit_length()} ==p:{w[u]==P}'
for u in (14853,1308,29967,11360,30163,24548,25442,7927,15616,11052,4432,19964):
    print(st(u))
print()
# who defines x_29967 / x_7927
for u in (29967, 7927):
    a=L.definer.get(u); print(f'definer of x_{u}: a{a} src={L.atom_src[a] if a is not None else None} neqs={len(L.atom2eq.get(a,{})) if a is not None else 0}')
print()
# --- repair a29539 after moving x_7068 ---
print('=== x_7068 move + a29539 repair ===')
for d in (1, P, 12345):
    w2=list(w); w2[7068]+=d; fwd2(w2)
    av=L.all_atom_values(w2)
    print(f' x_7068+{"p" if d==P else d}: a29539={av[29539]}   a29539 mod p={av[29539]%P}')
    changed={u:(w[u],w2[u]) for u in (14853,1308,29967,11360,30163) if w2[u]!=w[u]}
    print('   changed among a29539 vars/handles:', {k:(str(a)[:20],str(b)[:20]) for k,(a,b) in changed.items()})
    for h in (30163, 11360, 14853, 1308):
        w3=list(w2)
        nv=T.solve_lin(L.definer.get(29967) if False else None,0,w3) if False else None
        # solve a29539 = 0 via handle h by adjusting the defining atom of x_29967
        a2=L.definer.get(29967)
        # target x_29967 value:
        tgt = 12846437*(w3[14853]-w3[1308])
        # x_29967 = x_11360*x_30163  -> set h
        other = 11360 if h==30163 else (30163 if h==11360 else None)
        if other is None: continue
        if w3[other]==0: print(f'   handle x_{h}: other factor is 0, cannot'); continue
        if tgt % w3[other]: print(f'   handle x_{h}: NOT divisible (need {tgt%w3[other]!=0})'); continue
        w3[h]=tgt//w3[other]; fwd2(w3)
        av3=L.all_atom_values(w3); f3=set(L.failing_eqs(av3))
        print(f'   handle x_{h} -> score {L.NEQ-len(f3)}  out12={len(f3-E)} {sorted(f3-E)[:8]}')

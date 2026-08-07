"""Is a failing equation PERMANENTLY unfixable?
   If value(e) != 0 mod p and EVERY variable's delta on e is = 0 mod p, then e can never be
   satisfied from this state -- no enrichment of the move set helps."""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import atomval, load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
for path,nm in [(os.path.join(LAB,'best','new_instance_partial_39026.json'),'checkpoint 39026'),
                (os.path.join(HERE,'data','finish3_named.json'),'s11 best 39018')]:
    v=load_raw(path)
    AV=[atomval(a,v) for a in range(L.NA)]
    def eqs(e): return sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())
    FAIL=[e for e in range(L.NEQ) if eqs(e)!=0]
    print(f"=== {nm}: {len(FAIL)} failing")
    perm=0
    for e in FAIL:
        val=eqs(e)
        vs=set()
        for a in L.eq_atoms[e][2]: vs |= set(L.avars[a])
        live=[]
        for u in vs:
            old=v[u]; v[u]=old+1
            d=0
            for a in L.var_atoms[u]:
                if e in L.atom2eq.get(a,{}):
                    d += L.atom2eq[a][e]*(atomval(a,v)-AV[a])
            v[u]=old
            if d % P: live.append(u)
        tag = 'PERMANENT' if (val%P!=0 and not live) else ('fixable-modp' if live else 'value=0 modp')
        if val%P!=0 and not live: perm+=1
        print(f"   eq {e:6d}: val%p!=0 {val%P!=0}  vars={len(vs)}  vars with mod-p leverage={len(live)}  -> {tag}")
    print(f"   PERMANENTLY unfixable: {perm} of {len(FAIL)}")

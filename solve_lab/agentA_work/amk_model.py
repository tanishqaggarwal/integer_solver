"""Build the exact enlarged region model: atoms A*, knobs K, equations R*,
and the exact polynomial of each atom in the knobs (others frozen)."""
import sys, collections, json; sys.path.insert(0,'.')
import env, lib as L
P=env.P
v0=env.load_best(); av0=L.all_atom_values(v0)
E=sorted(set(e for a in env.SEVEN for e in L.atom2eq[a]))
A0=set(a for e in E for a in L.eq_atoms[e][2])

def build(extra):
    A=set(A0)|set(extra)
    K=sorted(u for u in set(u for a in A for u in L.avars[a]) if all(x in A for x in L.var_atoms[u]))
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    return sorted(A),K,R

def knobpoly(a,K,v):
    """polys[a] with non-knob vars frozen at v: dict monomial(tuple of knob idx)->coeff"""
    ki={u:i for i,u in enumerate(K)}
    out=collections.defaultdict(int)
    for m,c in L.polys[a].items():
        t=c; mono=[]
        for u in m:
            if u in ki: mono.append(ki[u])
            else: t*=v[u]
        if t: out[tuple(sorted(mono))]+=t
    return {k:c for k,c in out.items() if c}

if __name__=='__main__':
    for extra in ([],[37887],[37887,41906],[37887,41906,29426]):
        A,K,R=build(extra)
        print('extra=%-22s atoms=%d knobs=%d eqs=%d'%(str(extra),len(A),len(K),len(R)))
    A,K,R=build([37887,41906])
    print('knobs:',K)
    print()
    for a in A:
        kp=knobpoly(a,K,v0)
        deg=max((len(m) for m in kp),key=lambda x:x) if kp else 0
        deg=max([len(m) for m in kp]+[0])
        nl=[m for m in kp if len(m)>1]
        print('a%-6d  terms=%-4d degree_in_knobs=%d  nonlinear_terms=%d  const=%s'%(
            a,len(kp),deg,len(nl),'0' if kp.get((),0)==0 else 'NZ'))

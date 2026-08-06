"""Lift a mod-p solution to an exact integer assignment.

Once every atom is 0 mod p, each check's integer value is a multiple of p, and each check's
handle enters with coefficient exactly +-p (through a wire).  So the integer repair that failed
all session -- 'not divisible' -- must now succeed: that failure WAS the mod-p obstruction.

  1. lift every free input to its residue in [0,p)
  2. forward-evaluate over Z so all gates are exact
  3. for each failing check, solve it exactly for a free variable occurring linearly
  4. iterate to a fixpoint
"""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
import fw
P=L.P; sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)

def lin_solve(a,t,v):
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc
        for u in m:
            if u!=t: term*=v[u]
        c+=term
    if c==0: return None
    old=v[t]; v[t]=0; rest=L.evalpoly(L.polys[a],v); v[t]=old
    if rest%c: return None
    return -rest//c

def forward_frozen(v, frozen):
    for comp in fw.ORDER:
        if len(comp)==1:
            u=comp[0]
            if u in frozen: continue
            x=fw.solve_lin(L.definer[u],u,v)
            if x is not None: v[u]=x
        else:
            for _ in range(40):
                ch=False
                for u in comp:
                    if u in frozen: continue
                    x=fw.solve_lin(L.definer[u],u,v)
                    if x is not None and x!=v[u]: v[u]=x; ch=True
                if not ch: break
    return v

def report(v,tag):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]
    Bg=[a for a in B if L.atom_out.get(a) is not None]
    Bc=[a for a in B if L.atom_out.get(a) is None]
    F=L.failing_eqs(AV)
    print(f"{tag}: broken gates={len(Bg)} {Bg[:8]} broken checks={len(Bc)} {Bc[:8]} "
          f"failing eqs={len(F)} score={L.NEQ-len(F)}", flush=True)
    return Bc,F

def lift(modp, frozen=(), rounds=40):
    v=[int(x)%P for x in modp]
    frozen=set(frozen)
    forward_frozen(v, frozen)
    report(v,'  after Z forward eval')
    for it in range(rounds):
        AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
        bad=[a for a in range(L.NA) if AV[a]!=0 and L.atom_out.get(a) is None]
        if not bad: break
        did=False
        for a in bad:
            for t in sorted(L.avars[a], key=lambda t: len(L.var_atoms[t])):
                if t not in FREE or t in frozen: continue
                x=lin_solve(a,t,v)
                if x is None or x==v[t]: continue
                v[t]=x; forward_frozen(v,frozen); did=True; break
            if did: break
        if not did:
            print(f"  stuck at it{it}: {len(bad)} broken checks {bad[:10]}"); break
    return v

if __name__=='__main__':
    src=sys.argv[1]
    frozen=set(int(x) for x in sys.argv[2:])
    modp=json.load(open(src))
    print(f"lifting {os.path.basename(src)} (frozen {sorted(frozen)})")
    v=lift(modp, frozen)
    Bc,F=report(v,'  RESULT')
    out=os.path.join(HERE,'data','lift_out.json')
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(out,'w'))
    print("saved",out)

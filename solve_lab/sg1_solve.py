"""Damped coordinate-descent fail-minimizer with incremental eval. Monotone: accept only
fail-reducing (or fail-equal + residual-reducing) integer moves. Knobs = frees in failing eqs."""
import sys, json, time, random, math
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p; val = H.val
from collections import defaultdict
ns = {'v': val, '__builtins__': {}}
eqcode = H.eqcode; eqvars = H.eqvars

# free -> eqs it can affect (through gates)
efs = []
for i in range(len(eqvars)):
    s=set()
    for var in eqvars[i]:
        if var in H.freeinp: s.add(var)
        else: s |= H.anc.get(var,set())
    efs.append(s)
free_to_eqs = defaultdict(set)
for i,s in enumerate(efs):
    for f in s: free_to_eqs[f].add(i)

sm = json.load(open('sg1_slavemap.json'))
checked = set(int(v) for v in sm)

def load(path):
    d = H.loadd(path)
    for v in H.freeinp: val[v] = d.get(v,0)
    H.forward()

def all_fails():
    return set(i for i,c in enumerate(eqcode) if eval(c,ns)!=0)

def eq_val(i):
    return eval(eqcode[i], ns)

def total_residual(F):
    return sum(abs(eval(eqcode[i],ns)) for i in F)

def candidates_for(f, F):
    """quadratic-root candidate values of f that zero some failing eq (via 3-pt sampling w/ forward)."""
    cands=set()
    x0=val[f]
    samp={}
    for dx in (0,1,2,-1,-2):
        val[f]=x0+dx; H.forward()
        samp[dx]=[eval(eqcode[i],ns) for i in F]
    val[f]=x0; H.forward()
    Flist=list(F)
    for k,i in enumerate(Flist):
        y0=samp[0][k]; y1=samp[1][k]; y2=samp[2][k]
        # fit y = A dx^2 + B dx + C
        C=y0; A=(y2-2*y1+y0)//2 if (y2-2*y1+y0)%2==0 else None; B=y1-y0-(A if A else 0)
        if A is None: 
            # linear fallback
            if y1!=y0 and (y0)%(y1-y0)==0: cands.add(x0 - y0//(y1-y0))
            continue
        if A==0:
            if B!=0 and (-C)%B==0: cands.add(x0 + (-C)//B)
        else:
            disc=B*B-4*A*C
            if disc>=0:
                r=math.isqrt(disc)
                if r*r==disc:
                    for s in (r,-r):
                        num=-B+s
                        if num%(2*A)==0: cands.add(x0+num//(2*A))
    return cands

def greedy(F0, knobs, budget_s=60, verbose=True):
    F=set(F0); t0=time.time(); rounds=0
    while time.time()-t0 < budget_s:
        rounds+=1; improved=False
        random.shuffle(knobs)
        for f in knobs:
            relF = F & free_to_eqs[f]
            if not relF: continue
            cands = candidates_for(f, relF)
            x0=val[f]; bestx=x0; bestn=len(F); bestr=total_residual(F)
            for x in cands:
                if x==x0: continue
                val[f]=x; H.forward()
                # incremental: recompute only free_to_eqs[f]
                newF = set(F)
                for i in free_to_eqs[f]:
                    if eval(eqcode[i],ns)!=0: newF.add(i)
                    else: newF.discard(i)
                n=len(newF)
                if n<bestn or (n==bestn and total_residual(newF)<bestr):
                    bestn=n; bestx=x; bestr=total_residual(newF); bestF=newF
            val[f]=bestx; H.forward()
            if bestx!=x0 and bestn<=len(F):
                if bestn<len(F): 
                    F=bestF; improved=True
                    if verbose: print(f"  round {rounds}: f={f} -> {len(F)} fails", flush=True)
        if not improved:
            if verbose: print(f"  plateau at {len(F)} after {rounds} rounds", flush=True)
            break
    return F

if __name__ == '__main__':
    start = sys.argv[1] if len(sys.argv)>1 else 'best_agentA_39022.json'
    budget = int(sys.argv[2]) if len(sys.argv)>2 else 60
    load(start)
    F0=all_fails()
    print(f"start {start}: {len(F0)} fails: {sorted(F0)}")
    # knobs: unchecked frees in support of failing eqs
    knobs=set()
    for i in F0: knobs |= (efs[i] - checked)
    knobs=sorted(knobs)
    print(f"knobs (unchecked frees in fail support): {len(knobs)}: {knobs[:40]}")
    F=greedy(F0, knobs, budget_s=budget)
    print(f"FINAL: {len(F)} fails: {sorted(F)}")
    if len(F)<len(F0):
        json.dump({f"x_{i}":val[i] for i in range(H.NVARS)}, open('sg1_greedy_out.json','w'))
        print("saved sg1_greedy_out.json")

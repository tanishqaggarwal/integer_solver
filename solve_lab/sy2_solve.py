import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
fc=H.loadd('fc_partial.json')
COMPOSITE={8731:{8731:1,4432:1}, 9118:{9118:1,7068:1}}
FIXED={4287,2081,24601,31861,14865,4432,7068}
def setup():
    for v in H.freeinp: H.val[v]=fc.get(v,0)
    H.forward()
def slave():
    H.val[4432]=H.val[19964]+H.val[28730]
    H.val[7068]=7376877*H.val[642]+H.val[2099]
    H.forward()
setup(); slave()
def build_closure(F):
    clo=set(F)
    for _ in range(4):
        fr=set()
        for i in clo: fr|=(H.eqvars[i]&H.freeinp)
        fr-=FIXED
        new=set()
        for kn in fr:
            for i,vs in enumerate(H.eqvars):
                if kn in vs: new.add(i)
        if new<=clo: break
        clo|=new
        if len(clo)>1600: break
    return sorted(clo)
def resids(Feqs):
    ns={'v':H.val,'__builtins__':{}}
    return [eval(H.eqcode[i],ns)%p for i in Feqs]
def apply(mv,s=1):
    for k,v in mv.items(): H.val[k]=H.val[k]+s*v
def modp_newton(iters=8):
    for it in range(iters):
        F=H.fails()
        if not F: return True
        Feqs=build_closure(F)
        base=resids(Feqs)
        nz=[k for k in range(len(Feqs)) if base[k]!=0]
        if not nz: return True   # all ≡0 mod p (only integer-lift left)
        knobs=set()
        for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
        knobs-=FIXED; knobs=sorted(knobs)
        cols={}
        for kn in knobs:
            mv=COMPOSITE.get(kn,{kn:1})
            apply(mv,1); H.forward(); Rn=resids(Feqs); apply(mv,-1); H.forward()
            c={k:(Rn[k]-base[k])%p for k in nz if (Rn[k]-base[k])%p!=0}
            if c: cols[kn]=c
        knobs=[k for k in cols]
        rowdata={k:{} for k in nz}
        for kn in knobs:
            for k,v in cols[kn].items(): rowdata[k][kn]=v
        rhs={k:(-base[k])%p for k in nz}
        used=set(); delta={}
        for kn in knobs:
            prow=None
            for r in nz:
                if r in used: continue
                if rowdata[r].get(kn,0)%p!=0: prow=r;break
            if prow is None: continue
            used.add(prow)
            ipv=inv(rowdata[prow][kn])
            for c in list(rowdata[prow]): rowdata[prow][c]=(rowdata[prow][c]*ipv)%p
            rhs[prow]=(rhs[prow]*ipv)%p
            for r in nz:
                if r==prow: continue
                f=rowdata[r].get(kn,0)%p
                if f==0: continue
                for c,val in rowdata[prow].items(): rowdata[r][c]=(rowdata[r].get(c,0)-f*val)%p
                rhs[r]=(rhs[r]-f*rhs[prow])%p
        incon=[r for r in nz if r not in used and rhs[r]%p!=0 and all(v%p==0 for v in rowdata[r].values())]
        # back-substitution: delta[kn] = rhs[pivrow] - sum over non-pivot cols... 
        # since we normalized+eliminated, pivot rows now: delta[pivkn] + sum(nonpiv coeffs*delta_nonpiv)=rhs
        # set non-pivot deltas=0 -> delta[pivkn]=rhs[pivrow]
        pivkn={}
        for kn in knobs:
            for r in nz:
                if rowdata[r].get(kn,0)%p==1 and r in used:
                    # check it's the pivot (only this row normalized to 1 for kn)
                    pass
        # simpler: recompute pivot mapping
        used2=set(); delta={}
        # redo elimination bookkeeping cleanly by tracking pivot per column
        # (re-run minimal): assign delta from pivot rows
        # Since we eliminated in-place, each pivot row r has rowdata[r][kn]=1 for its pivot kn and 0 for other pivots
        # Identify pivot kn for each used row:
        for r in used:
            piv=[kn for kn in knobs if rowdata[r].get(kn,0)%p!=0]
            # pivot col is the one unique to this row
        # fallback: solve by taking delta[kn]=rhs[prow] using stored order
        # -- rebuild via explicit pivot tracking
        return ('need_explicit',incon,Feqs)
    return False
# The above got messy; do explicit pivot tracking version:
def modp_newton2(iters=10):
    for it in range(iters):
        F=H.fails()
        if not F: 
            print('  iter',it,'ALL PASS'); return True
        Feqs=build_closure(F)
        base=resids(Feqs)
        nz=[k for k in range(len(Feqs)) if base[k]!=0]
        if not nz:
            print('  iter',it,'residuals all ≡0 mod p (lift remains)'); return 'lift'
        knobs=set()
        for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
        knobs-=FIXED; knobs=sorted(knobs)
        cols={}
        for kn in knobs:
            mv=COMPOSITE.get(kn,{kn:1})
            apply(mv,1); H.forward(); Rn=resids(Feqs); apply(mv,-1); H.forward()
            c={k:(Rn[k]-base[k])%p for k in nz if (Rn[k]-base[k])%p!=0}
            if c: cols[kn]=c
        knobs=[k for k in cols]
        rowdata={k:dict((kn,cols[kn][k]) for kn in knobs if k in cols[kn]) for k in nz}
        rhs={k:(-base[k])%p for k in nz}
        pivot=[]  # (kn,row)
        usedrows=set()
        for kn in knobs:
            prow=None
            for r in nz:
                if r in usedrows: continue
                if rowdata[r].get(kn,0)%p!=0: prow=r;break
            if prow is None: continue
            usedrows.add(prow); pivot.append((kn,prow))
            ipv=inv(rowdata[prow][kn])
            for c in list(rowdata[prow]): rowdata[prow][c]=(rowdata[prow][c]*ipv)%p
            rhs[prow]=(rhs[prow]*ipv)%p
            for r in nz:
                if r==prow: continue
                f=rowdata[r].get(kn,0)%p
                if f==0: continue
                for c,val in rowdata[prow].items(): rowdata[r][c]=(rowdata[r].get(c,0)-f*val)%p
                rhs[r]=(rhs[r]-f*rhs[prow])%p
        incon=[r for r in nz if r not in usedrows and rhs[r]%p!=0]
        # solution: non-pivot delta=0 -> delta[pivkn]=rhs[pivrow]
        delta={kn:rhs[prow] for kn,prow in pivot}
        for kn,dv in delta.items():
            mv=COMPOSITE.get(kn,{kn:1})
            apply({k:v*dv for k,v in mv.items()},1)
        H.forward(); slave()
        nf=len(H.fails())
        print('  iter %d: rows=%d rank=%d incon=%d -> fails=%d'%(it,len(nz),len(pivot),len(incon),nf))
    return False
print('start fails:',len(H.fails()))
res=modp_newton2()
print('modp result:',res,' fails now:',len(H.fails()))
# integer lift: set p-granular handles for the 6 gadgets
def lift():
    H.forward()
    if H.val[9106]%(13523997*p)==0: H.val[950]=H.val[9106]//(13523997*p)
    if H.val[2239]%p==0: H.val[6947]=(6122989*H.val[2239])//p
    if H.val[31731]%p==0: H.val[33168]=-(H.val[31731]//p)
    H.forward(); slave()
lift()
F=H.fails()
print('after lift: fails=',len(F), sorted(F)[:20])
# save
out={f'x_{i}':H.val[i] for i in range(H.NVARS) if H.val[i]!=0}
json.dump(out,open('sy2_attempt.json','w'))

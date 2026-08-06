import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
# Build config: regime (1,1) with load-condition residues
r8=109044024666698959972204451600908701898659086097062528124234304603594878834481
r9=33371159155735472537534252650716501592825364489306217536352743247010353604716
vA=H.loadd('best_agentA_39022.json')
def setup():
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    H.val[4287]=1
    H.val[31861]=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
    H.val[14865]=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
    H.val[8731]=r8; H.val[9118]=r9
    H.val[9413]=0; H.val[17325]=0
    H.forward()
    H.val[4432]=H.val[19964]; H.val[7068]=H.val[2099]
    H.forward()
    # load handles
    H.val[950]=H.val[9106]//(13523997*p)
    H.val[6947]=(6122989*H.val[2239])//p
    H.val[33168]=-(H.val[31731]//p)
    H.forward()
setup()
F=H.fails()
print('start fails:',len(F))
# relevant free inputs = free inputs in the union of failing eqs
def relevant(F):
    s=set()
    for i in F: s|=(H.eqvars[i]&H.freeinp)
    return sorted(s)
# mod-p Newton: linearize residuals of F wrt knobs, solve delta mod p
def residsF(F):
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns)%p for i in F}
def mod_newton(iters=6):
    for it in range(iters):
        F=H.fails()
        if not F:
            print('  ALL PASS'); return
        R=residsF(F)
        Fp=[i for i in F if R[i]!=0]      # nonzero mod p
        # some fails may be ≡0 mod p (pure p-lift) - skip in mod phase
        knobs=relevant(F)
        # build jacobian columns
        base={i:R[i] for i in F}
        cols={}
        for kn in knobs:
            old=H.val[kn]; H.val[kn]=(old+1); H.forward()
            Rn=residsF(F)
            H.val[kn]=old; H.forward()
            col={i:(Rn[i]-base[i])%p for i in F if (Rn[i]-base[i])%p!=0}
            if col: cols[kn]=col
        # Gaussian elimination over GF(p): solve sum_kn delta[kn]*col = -base
        # Represent as rows=eqs, vars=knobs
        knobs=[k for k in cols]
        rowset=sorted(set(i for c in cols.values() for i in c) | set(base))
        # build dense matrix
        A={(i,kn):cols[kn].get(i,0) for kn in knobs for i in rowset}
        bvec={i:(-base[i])%p for i in rowset}
        # gaussian
        delta={kn:0 for kn in knobs}
        rows=rowset[:]; used=set()
        pivcols=[]
        Aw={i:{kn:A[(i,kn)]%p for kn in knobs} for i in rows}
        bw=dict(bvec)
        colorder=knobs[:]
        pivot_for={}
        r=0
        rlist=rows[:]
        for kn in colorder:
            # find pivot row
            prow=None
            for i in rlist:
                if i in used: continue
                if Aw[i].get(kn,0)%p!=0: prow=i;break
            if prow is None: continue
            used.add(prow); pivot_for[kn]=prow
            iv=inv(Aw[prow][kn])
            for i in rlist:
                if i==prow or i in used and i!=prow: 
                    pass
            # eliminate kn from all other rows
            for i in rlist:
                if i==prow: continue
                f=Aw[i].get(kn,0)%p
                if f==0: continue
                for kn2 in colorder:
                    Aw[i][kn2]=(Aw[i][kn2]-f*iv*Aw[prow].get(kn2,0))%p
                bw[i]=(bw[i]-f*iv*bw[prow])%p
        # back-substitute (diagonal-ish): for each pivot col, delta = bw[prow]*inv(piv)
        for kn,prow in pivot_for.items():
            delta[kn]=(bw[prow]*inv(Aw[prow][kn]))%p
        # check consistency: rows with all-zero A but nonzero b => inconsistent
        incon=[i for i in rlist if i not in used and bw[i]%p!=0 and all(Aw[i].get(kn,0)%p==0 for kn in colorder)]
        # apply
        for kn,dv in delta.items():
            H.val[kn]=(H.val[kn]+dv)%p
        H.forward()
        nf=len(H.fails())
        print(f'  iter {it}: knobs={len(knobs)} pivots={len(pivot_for)} inconsistent_rows={len(incon)} -> fails={nf}')
        if incon: 
            print('    INCONSISTENT mod-p rows:',incon[:8]); 
    return
mod_newton()
F=H.fails()
print('after mod-newton: fails=',len(F), sorted(F)[:30])

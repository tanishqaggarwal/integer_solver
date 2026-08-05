import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
fc=H.loadd('fc_partial.json')
COMPOSITE={8731:{8731:1,4432:1}, 9118:{9118:1,7068:1}}
FIXED={4287,2081,24601,31861,14865,4432,7068}
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.forward()
def slave():
    H.val[4432]=H.val[19964]+H.val[28730]; H.val[7068]=7376877*H.val[642]+H.val[2099]; H.forward()
slave()
# closure to convergence
clo=set(H.fails())
while True:
    fr=set()
    for i in clo: fr|=(H.eqvars[i]&H.freeinp)
    fr-=FIXED
    new=set()
    for kn in fr:
        for i,vs in enumerate(H.eqvars):
            if kn in vs: new.add(i)
    if new<=clo or len(clo)>2500: break
    clo|=new
Feqs=sorted(clo)
knobs=set()
for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
knobs-=FIXED; knobs=sorted(knobs)
print('closure eqs:',len(Feqs),' knobs:',len(knobs))
# eq index map
pos={e:k for k,e in enumerate(Feqs)}
# precompute which closure eqs each knob (composite parts) touches
knob_eqs={}
for kn in knobs:
    parts=set(COMPOSITE.get(kn,{kn:1}))
    knob_eqs[kn]=[i for i in Feqs if H.eqvars[i]&parts]
def eval_eqs(idxs):
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns)%p for i in idxs}
def newton(iters=8):
    for it in range(iters):
        F=H.fails()
        allr=eval_eqs(Feqs)
        nz=[i for i in Feqs if allr[i]!=0]
        if not nz:
            print('  iter',it,': all closure residuals ≡0 mod p'); return 'lift'
        # build sparse columns (mod p)
        cols={}
        for kn in knobs:
            mv=COMPOSITE.get(kn,{kn:1})
            for k,v in mv.items(): H.val[k]+=v
            H.forward()
            ns={'v':H.val,'__builtins__':{}}
            col={}
            for i in knob_eqs[kn]:
                d=(eval(H.eqcode[i],ns)-allr[i])%p
                if d: col[i]=d
            for k,v in mv.items(): H.val[k]-=v
            H.forward()
            if col: cols[kn]=col
        # sparse GF(p) gaussian over ALL Feqs rows
        rowdata={i:{} for i in Feqs}
        for kn,col in cols.items():
            for i,v in col.items(): rowdata[i][kn]=v
        rhs={i:(-allr[i])%p for i in Feqs}
        kl=[k for k in knobs if k in cols]
        used=set(); pivots=[]
        for kn in kl:
            prow=None
            for i in Feqs:
                if i in used: continue
                if rowdata[i].get(kn,0)%p!=0: prow=i;break
            if prow is None: continue
            used.add(prow); pivots.append((kn,prow))
            ipv=inv(rowdata[prow][kn])
            for c in list(rowdata[prow]): rowdata[prow][c]=(rowdata[prow][c]*ipv)%p
            rhs[prow]=(rhs[prow]*ipv)%p
            for i in Feqs:
                if i==prow: continue
                f=rowdata[i].get(kn,0)%p
                if f==0: continue
                for c,val in rowdata[prow].items(): rowdata[i][c]=(rowdata[i].get(c,0)-f*val)%p
                rhs[i]=(rhs[i]-f*rhs[prow])%p
        incon=[i for i in Feqs if i not in used and rhs[i]%p!=0]
        delta={kn:rhs[prow] for kn,prow in pivots}
        print('  iter %d: nz=%d rank=%d incon=%d'%(it,len(nz),len(pivots),len(incon)))
        if incon:
            print('    INCONSISTENT eqs:',incon[:12]); return ('incon',incon)
        for kn,dv in delta.items():
            mv=COMPOSITE.get(kn,{kn:1})
            for k,v in mv.items(): H.val[k]+=v*dv
        H.forward(); slave()
    return 'maxiter'
r=newton()
print('newton:',r,' fails now:',len(H.fails()))
# lift
H.forward()
if H.val[9106]%(13523997*p)==0: H.val[950]=H.val[9106]//(13523997*p)
if H.val[2239]%p==0: H.val[6947]=(6122989*H.val[2239])//p
if H.val[31731]%p==0: H.val[33168]=-(H.val[31731]//p)
H.forward(); slave()
F=H.fails()
print('after lift exact fails:',len(F), sorted(F)[:20])
if len(F)<11:
    out={f'x_{i}':H.val[i] for i in range(H.NVARS) if H.val[i]!=0}
    json.dump(out,open('sy2_best.json','w')); print('SAVED sy2_best.json',len(F),'fails')

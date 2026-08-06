import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
fc=H.loadd('fc_partial.json')
PROTECT={4287,2081,24601,31861,14865,6418,12553,8731,9118,4432,7068}
def setup():
    for v in H.freeinp: H.val[v]=fc.get(v,0)
    H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0; H.val[950]=0; H.val[6947]=0
    H.forward()
    if H.val[37720]%9994531==0: H.val[8976]=H.val[37720]//9994531
    H.forward()
setup()
print('start fails:',len(H.fails()),sorted(H.fails()))
# closure from consumer fails
clo=set(H.fails())
while True:
    fr=set()
    for i in clo: fr|=(H.eqvars[i]&H.freeinp)
    fr-=PROTECT
    new=set()
    for kn in fr:
        for i,vs in enumerate(H.eqvars):
            if kn in vs: new.add(i)
    if new<=clo or len(clo)>2000: break
    clo|=new
Feqs=sorted(clo)
knobs=set()
for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
knobs-=PROTECT; knobs=sorted(knobs)
print('closure eqs:',len(Feqs),' knobs:',len(knobs))
knob_eqs={kn:[i for i in Feqs if kn in H.eqvars[i]] for kn in knobs}
def eval_all():
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns) for i in Feqs}
def newton(iters=10):
    for it in range(iters):
        F=H.fails()
        if not F: print('  ALL PASS at iter',it); return True
        allr=eval_all()
        base={i:allr[i]%p for i in Feqs}
        nz=[i for i in Feqs if base[i]!=0]
        if not nz: print('  iter',it,'all ≡0 mod p'); return 'lift'
        cols={}
        for kn in knobs:
            H.val[kn]+=1; H.forward()
            ns={'v':H.val,'__builtins__':{}}
            col={}
            for i in knob_eqs[kn]:
                d=(eval(H.eqcode[i],ns)%p - base[i])%p
                if d: col[i]=d
            H.val[kn]-=1; H.forward()
            if col: cols[kn]=col
        rowdata={i:{} for i in Feqs}
        for kn,col in cols.items():
            for i,v in col.items(): rowdata[i][kn]=v
        rhs={i:(-base[i])%p for i in Feqs}
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
        print('  iter %d: nz=%d rank=%d incon=%d'%(it,len(nz),len(pivots),len(incon)))
        if incon: print('    INCON:',incon[:10]); return ('incon',incon)
        for kn,prow in pivots:
            H.val[kn]+=rhs[prow]
        H.forward()
    return 'maxiter'
r=newton()
print('newton:',r,'fails now:',len(H.fails()))
# integer lift: p-granular slacks. Iterate: for each nonzero eq, find p-granular unit-of-p slack
F=H.fails()
print('after modp, exact fails:',len(F),sorted(F)[:20])
if len(F)<11:
    out={f'x_{i}':H.val[i] for i in range(H.NVARS) if H.val[i]!=0}
    json.dump(out,open('sy2_heal_out.json','w')); print('SAVED',len(F))

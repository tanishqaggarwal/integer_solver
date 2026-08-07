"""Per-node coordinate alignment (from slot links) + per-node chord orientation (numeric)."""
import sys, os, json, collections, pickle, re, time
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
D=pickle.load(open('ortree2.pkl','rb')); tree=dict(D['tree'])
NODE=dict(pickle.load(open('nodes.pkl','rb'))); OUT=dict(pickle.load(open('outwires.pkl','rb')))
P=pickle.load(open('pins.pkl','rb')); PIN=P['PIN']; live=P['live']; dead=P['dead']
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
def deref(v):
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen: seen.add(v); v=defrhs[v][1]
    return v
RA,RB=D['RA'],D['RB']
prodof=collections.defaultdict(list); notof=collections.defaultdict(list); gated=collections.defaultdict(list); sumof={}
for w,r in defrhs.items():
    if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v':
        a,b=deref(r[1][1]),deref(r[2][1]); prodof[frozenset((a,b))].append(w)
        gated[a].append((w,b)); gated[b].append((w,a))
    elif r[0]=='-' and r[1][0]=='c' and r[1][1]==1 and r[2][0]=='v': notof[deref(r[2][1])].append(deref(w))
    elif r[0]=='+' and r[1][0]=='v' and r[2][0]=='v': sumof[frozenset((r[1][1],r[2][1]))]=w
ROOT=dict(a=RA,b=RB,sa=[x for nb in notof[RB] for x in prodof.get(frozenset((RA,nb)),[])],
          sb=[x for na in notof[RA] for x in prodof.get(frozenset((RB,na)),[])],
          sab=prodof.get(frozenset((RA,RB)),[]))
for k,s in (('ga','sa'),('gb','sb'),('gab','sab')): ROOT[k]=[g for x in ROOT[s] for g in gated[x]]
def outwire(N):
    res=[]
    for (pa,va) in N['ga']:
        c=[(pb,vb) for (pb,vb) in N['gb'] if frozenset((pa,pb)) in sumof]
        if len(c)!=1: return None
        pb,vb=c[0]; s2=sumof[frozenset((pa,pb))]
        c2=[(pab,vab) for (pab,vab) in N['gab'] if frozenset((s2,pab)) in sumof]
        if len(c2)!=1: return None
        pab,vab=c2[0]; res.append(dict(va=va,vb=vb,vab=vab,out=sumof[frozenset((s2,pab))]))
    return res if len(res)==2 else None
ROOTID=-1
NODE[ROOTID]=ROOT; OUT[ROOTID]=outwire(ROOT); tree[ROOTID]=(RA,RB)
link={}
lr=[re.compile(r'^\(\(x(\d+)-x(\d+)\)[-+]'), re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-')]
for a in E.res:
    m=lr[0].match(a)
    if m: u,z=int(m.group(1)),int(m.group(2))
    else:
        m=lr[1].match(a)
        if not m: continue
        u,z=int(m.group(2)),int(m.group(3))
    if u not in defrhs and z in defrhs: link[u]=z
    elif z not in defrhs and u in defrhs: link[z]=u
# ---- coordinate alignment ----
perm={}   # (n,side) -> [j0,j1] : parent coord i corresponds to child coord perm[i]
bad=[]
for n in NODE:
    a,b=tree[n]
    for side,ch in (('va',a),('vb',b)):
        if tree[ch] is None: perm[(n,side)]=[0,1]; continue
        couts=[d['out'] for d in OUT[ch]]
        pm=[]
        for i,d in enumerate(OUT[n]):
            z=link.get(d[side])
            pm.append(couts.index(z) if z in couts else None)
        if None in pm or sorted(pm)!=[0,1]: bad.append((n,side,pm))
        perm[(n,side)]=pm

# ---- numeric repair of unresolved coordinate alignments ----
H=pickle.load(open('handles.pkl','rb')); appearP=H['appearP']
print('coordinate alignment: %d node-sides, %d bad (pre-repair)'%(len(perm),len(bad)))
for x in bad[:5]: print('   ',x)
# ---- subtree leaf sets ----
sub={}
def subl(n):
    if n in sub: return sub[n]
    if tree[n] is None: sub[n]=[n]
    else: sub[n]=subl(tree[n][0])+subl(tree[n][1])
    return sub[n]
subl(ROOTID)
liveset=set(live)
# ---- chord orientation, numerically, per node ----
def chord(A,B,orient):
    ax,ay,bx,by=(A[orient],A[1-orient],B[orient],B[1-orient])
    if (bx-ax)%p==0: return None
    l=(by-ay)*pow(bx-ax,p-2,p)%p
    ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    return (ox,oy) if orient==0 else (oy,ox)
ORIENT={}; VALCACHE={}
def buildvals(S, ORI):
    """return (isl,valn) with valn[n] in node n's own coord order"""
    isl={}; valn={}
    order=[]
    def post(n):
        if tree[n] is None:
            isl[n]= n in S; valn[n]= tuple(PIN[n][1]) if n in S else None; return
        for c in tree[n]: post(c)
        order.append(n)
    post(ROOTID)
    for n in order:
        a,b=tree[n]; la,lb=isl[a],isl[b]; isl[n]=la or lb
        def proj(ch,side):
            pm=perm[(n,side)]
            if valn[ch] is None or pm[0] is None or pm[1] is None: return None
            return (valn[ch][pm[0]],valn[ch][pm[1]])
        if la and lb:
            o=ORI.get(n)
            PA=proj(a,'va'); PB=proj(b,'vb')
            valn[n]=None if (o is None or PA is None or PB is None) else chord(PA,PB,o)
        elif la: valn[n]=proj(a,'va')
        elif lb: valn[n]=proj(b,'vb')
        else: valn[n]=None
    return isl,valn,order
def assignment(S,ORI):
    isl,valn,order=buildvals(S,ORI)
    v={}
    for L in S:
        ws,Cs=PIN[L]; v[L]=1
        for w,c in zip(ws,Cs): v[w]=c
    for n in order:
        a,b=tree[n]
        for i,d in enumerate(OUT[n]):
            for side,ch in (('va',a),('vb',b)):
                if tree[ch] is not None:
                    pmi=perm[(n,side)][i]
                    v[d[side]] = (valn[ch][pmi] if (isl[ch] and pmi is not None and valn[ch] is not None) else 0)
            v[d['vab']] = (valn[n][i] if (isl[a] and isl[b] and valn[n] is not None) else 0)
    return v,isl,valn,order
def run(v):
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    return vv,E.run(vv)

def repair_perm():
    fixed=0
    for (n,side),pm in list(perm.items()):
        if None not in pm: continue
        ch=tree[n][0 if side=='va' else 1]
        cand=[x for x in sub[ch] if x in set(live)]
        if not cand: continue
        S={cand[0]}
        v,isl,valn,_=assignment(S,{})
        target=valn[ch]
        newpm=list(pm)
        for i,d in enumerate(OUT[n]):
            if newpm[i] is not None: continue
            w=d[side]
            vv=dict(v); vv[w]=0; _,r0=run(vv)
            vv[w]=1; _,r1=run(vv)
            got=None
            for a in appearP.get(w,[]):
                idx=E.residx[a]; f0=r0[idx]%p; sl=(r1[idx]-r0[idx])%p
                if f0 and sl:
                    c=(-f0)*pow(sl,p-2,p)%p
                    if got is None: got=c
                    elif got!=c: got='CONFLICT'
            if got in target: newpm[i]=list(target).index(got)
        if None not in newpm and sorted(newpm)==[0,1]:
            perm[(n,side)]=newpm; fixed+=1
        else: print('   UNREPAIRED',n,side,newpm)
    print('perm repaired numerically:',fixed)
repair_perm()

if __name__=='__main__':
    t0=time.time()
    order=[]
    def post(n):
        if tree[n] is None: return
        for c in tree[n]: post(c)
        order.append(n)
    post(ROOTID)
    fails=[]
    for n in order:
        a,b=tree[n]
        la=[x for x in sub[a] if x in liveset]; lb=[x for x in sub[b] if x in liveset]
        if not la or not lb: ORIENT[n]='DEAD'; continue
        S={la[0],lb[0]}
        sel=NODE[n]['sab'][0]
        chk=[x for x in resby.get(sel,[]) if ('x%d*'%sel) in x or ('*x%d'%sel) in x]
        if len(chk)<3: fails.append((n,'checks',len(chk))); ORIENT[n]=None; continue
        ci=[E.residx[x] for x in chk[:3]]
        vabw=[d['vab'] for d in OUT[n]]
        v,isl,valn,_=assignment(S,ORIENT)
        v[vabw[0]]=0; v[vabw[1]]=0; _,r0=run(v); f0=[r0[i]%p for i in ci]
        v[vabw[0]]=1; _,r1=run(v); c1=[(r1[i]-r0[i])%p for i in ci]; v[vabw[0]]=0
        v[vabw[1]]=1; _,r2=run(v); c2=[(r2[i]-r0[i])%p for i in ci]; v[vabw[1]]=0
        det=(c1[0]*c2[1]-c1[1]*c2[0])%p
        if det==0: fails.append((n,'det0')); ORIENT[n]=None; continue
        di=pow(det,p-2,p)
        d0=((-f0[0])*c2[1]+f0[1]*c2[0])%p*di%p
        d1=(c1[0]*(-f0[1])+c1[1]*f0[0])%p*di%p
        third=(c1[2]*d0+c2[2]*d1+f0[2])%p==0
        A=(valn[a][perm[(n,'va')][0]],valn[a][perm[(n,'va')][1]])
        B=(valn[b][perm[(n,'vb')][0]],valn[b][perm[(n,'vb')][1]])
        got=(d0,d1); ok=None
        for o in (0,1):
            if chord(A,B,o)==got: ok=o
        if ok is None or not third: fails.append((n,'law',third,got,[chord(A,B,o) for o in (0,1)]))
        ORIENT[n]=ok
    print('nodes calibrated %d ; orientation hist %s ; failures %d'%(
        len(order),collections.Counter(ORIENT.values()).most_common(),len(fails)))
    for x in fails[:6]: print('   FAIL',x)
    print('elapsed %.0fs'%(time.time()-t0))
    pickle.dump({'ORIENT':ORIENT,'perm':perm,'tree':tree,'NODE':NODE,'OUT':OUT,'link':link,
                 'PIN':PIN,'live':live,'dead':dead,'sub':sub,'ROOTID':ROOTID,'order':order},
                open('calib.pkl','wb'))

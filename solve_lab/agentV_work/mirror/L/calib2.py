import sys, os, collections, pickle, time
F='/home/user/integer_solver/solve_lab/agentV_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
M=pickle.load(open('full_model.pkl','rb'))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; live=set(M['live']); link=M['link']
sub=M['sub']; order=M['order']; PIN=M['PIN']; ROOT=M['ROOT']
H=pickle.load(open('handles.pkl','rb')); appearP=H['appearP']
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
def run(v):
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    return vv,E.run(vv)
perm={}; bad=[]
for n in NODE:
    for side,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])):
        if tree[ch] is None: perm[(n,side)]=[0,1]; continue
        couts=[d['out'] for d in OUT[ch]]
        pm=[couts.index(link[d[side]]) if link.get(d[side]) in couts else None for d in OUT[n]]
        perm[(n,side)]=pm
        if None in pm or sorted(pm)!=[0,1]: bad.append((n,side,pm))
print('alignment bad (pre-repair):',len(bad))
def chord(A,B,o):
    ax,ay,bx,by=A[o],A[1-o],B[o],B[1-o]
    if (bx-ax)%p==0: return None
    l=(by-ay)*pow(bx-ax,p-2,p)%p
    ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    return (ox,oy) if o==0 else (oy,ox)
def buildvals(S,ORI):
    isl={}; valn={}
    for L in tree:
        if tree[L] is None: isl[L]=L in S; valn[L]=tuple(PIN[L][1]) if L in S else None
    for n in order:
        a,b=tree[n]; la,lb=isl[a],isl[b]; isl[n]=la or lb
        def proj(ch,side):
            pm=perm[(n,side)]
            if valn[ch] is None or pm[0] is None or pm[1] is None: return None
            return (valn[ch][pm[0]],valn[ch][pm[1]])
        if la and lb:
            o=ORI.get(n); PA=proj(a,'va'); PB=proj(b,'vb')
            valn[n]=None if (o is None or o=='DEAD' or PA is None or PB is None) else chord(PA,PB,o)
        elif la: valn[n]=proj(a,'va')
        elif lb: valn[n]=proj(b,'vb')
        else: valn[n]=None
    return isl,valn
def assignment(S,ORI):
    isl,valn=buildvals(S,ORI); v={}
    for L in S:
        ws,Cs=PIN[L]; v[L]=1
        for w,c in zip(ws,Cs): v[w]=c
    for n in order:
        a,b=tree[n]
        for i,d in enumerate(OUT[n]):
            for side,ch in (('va',a),('vb',b)):
                if tree[ch] is not None:
                    pmi=perm[(n,side)][i]
                    v[d[side]]=(valn[ch][pmi] if (isl[ch] and pmi is not None and valn[ch] is not None) else 0)
            v[d['vab']]=(valn[n][i] if (isl[a] and isl[b] and valn[n] is not None) else 0)
    return v,isl,valn
# numeric perm repair
for (n,side),pm in list(perm.items()):
    if None not in pm: continue
    ch=tree[n][0 if side=='va' else 1]
    cand=[x for x in sub[ch] if x in live]
    if not cand: continue
    v,isl,valn=assignment({cand[0]},{})
    tgt=valn[ch]; new=list(pm)
    for i,d in enumerate(OUT[n]):
        if new[i] is not None: continue
        w=d[side]; vv=dict(v); vv[w]=0; _,r0=run(vv); vv[w]=1; _,r1=run(vv)
        got=None
        for a in appearP.get(w,[]):
            j=E.residx[a]; f0=r0[j]%p; sl=(r1[j]-r0[j])%p
            if f0 and sl:
                c=(-f0)*pow(sl,p-2,p)%p
                got=c if got is None else ('X' if got!=c else got)
        if tgt and got in tgt: new[i]=list(tgt).index(got)
    if None not in new and sorted(new)==[0,1]: perm[(n,side)]=new
    else: print('  UNREPAIRED',n,side,new)
print('alignment bad (post-repair):',sum(1 for pm in perm.values() if None in pm or sorted(pm)!=[0,1]))
# orientation calibration
t0=time.time(); ORIENT={}; fails=[]
for n in order:
    a,b=tree[n]
    la=[x for x in sub[a] if x in live]; lb=[x for x in sub[b] if x in live]
    if not la or not lb: ORIENT[n]='DEAD'; continue
    sel=NODE[n]['sab'][0]
    chk=[x for x in resby.get(sel,[]) if ('x%d*'%sel) in x or ('*x%d'%sel) in x]
    if len(chk)<3: fails.append((n,'nchecks',len(chk))); ORIENT[n]=None; continue
    ci=[E.residx[x] for x in chk[:3]]; vw=[d['vab'] for d in OUT[n]]
    v,isl,valn=assignment({la[0],lb[0]},ORIENT)
    v[vw[0]]=0; v[vw[1]]=0; _,r0=run(v); f0=[r0[i]%p for i in ci]
    v[vw[0]]=1; _,r1=run(v); c1=[(r1[i]-r0[i])%p for i in ci]; v[vw[0]]=0
    v[vw[1]]=1; _,r2=run(v); c2=[(r2[i]-r0[i])%p for i in ci]; v[vw[1]]=0
    det=(c1[0]*c2[1]-c1[1]*c2[0])%p
    if det==0: fails.append((n,'det0')); ORIENT[n]=None; continue
    di=pow(det,p-2,p)
    d0=((-f0[0])*c2[1]+f0[1]*c2[0])%p*di%p
    d1=(c1[0]*(-f0[1])+c1[1]*f0[0])%p*di%p
    third=(c1[2]*d0+c2[2]*d1+f0[2])%p==0
    A=(valn[a][perm[(n,'va')][0]],valn[a][perm[(n,'va')][1]])
    B=(valn[b][perm[(n,'vb')][0]],valn[b][perm[(n,'vb')][1]])
    ok=None
    for o in (0,1):
        if chord(A,B,o)==(d0,d1): ok=o
    if ok is None or not third: fails.append((n,'law',third))
    ORIENT[n]=ok
print('calibrated %d nodes; orient hist %s; failures %d %s (%.0fs)'%(
    len(order),collections.Counter(ORIENT.values()).most_common(),len(fails),fails[:4],time.time()-t0))
pickle.dump({'ORIENT':ORIENT,'perm':perm},open('calib2.pkl','wb'))

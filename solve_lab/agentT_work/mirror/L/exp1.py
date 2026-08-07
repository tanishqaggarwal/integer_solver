"""EXPERIMENT 1 (agent L): two leaves ON in the SAME OR-group.
Measure by exact re-propagation what the slot wire actually holds."""
import sys, os, json, collections, pickle
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
C={10424:108269144428078338534087367253036776177727590608090039786410378426385300938364061820017016,
   27436:16230052175022190183706799464018524297465331898484712082891903119492101192039356374894798,
   20930:103584408103033301326011545465017242869465633026144536135655992693846868876792606063391405,
   30632:55691485823142661220390052154000291499156499018895928723006063316406569613942245404099691}
CHK=['((x17100*x35069)+x16531)','((9553893*(x17100*x18439))-x20490)','((x17100*x14819)-(7038713*x20150))']
ci=[E.residx[a] for a in CHK]
OUT=[35346,10824]                     # node x23153 mux output wires (coord1, coord2)
SELS={'sel_a':20753,'sel_b':2370,'sel_ab':17100}
GATED={'coord1':(10424,20930,36193),'coord2':(27436,30632,35256)}

def run(bits, vab=(0,0), pins=True):
    v=[0]*NV
    for l,b in bits.items(): v[l]=b
    if pins:
        for l,ws in ((19326,(10424,27436)),(28825,(20930,30632))):
            if bits.get(l): 
                for w in ws: v[w]=C[w]
    v[36193],v[35256]=vab
    r=E.run(v)
    return v,r

def solve_vab(bits):
    """solve the 3 stage checks for (x36193,x35256) mod p by affine probing"""
    _,r0=run(bits,(0,0)); f0=[r0[i]%p for i in ci]
    _,r1=run(bits,(1,0)); c1=[(r1[i]-r0[i])%p for i in ci]
    _,r2=run(bits,(0,1)); c2=[(r2[i]-r0[i])%p for i in ci]
    det=(c1[0]*c2[1]-c1[1]*c2[0])%p
    if det==0: return None,f0,c1,c2
    di=pow(det,p-2,p)
    d0=((-f0[0])*c2[1]+f0[1]*c2[0])%p*di%p
    d1=(c1[0]*(-f0[1])+c1[1]*f0[0])%p*di%p
    ok=(c1[2]*d0+c2[2]*d1+f0[2])%p==0
    return (d0,d1,ok),f0,c1,c2

def chordK(A,B):
    l=(B[1]-A[1])*pow(B[0]-A[0],p-2,p)%p
    ox=(l*l-A[0]-B[0]-K)%p
    oy=(l*(A[0]-ox)-A[1])%p
    return ox,oy

print('='*90)
for bits,label in [({19326:1,28825:0},'A only'),({19326:0,28825:1},'B only'),({19326:1,28825:1},'BOTH (same OR-group)')]:
    sol,f0,c1,c2=solve_vab(bits)
    print('--- %s ---'%label)
    v,r=run(bits,(sol[0],sol[1]) if sol else (0,0))
    print('  selectors: sel_a=%d sel_b=%d sel_ab=%d'%(v[20753]%p,v[2370]%p,v[17100]%p))
    print('  chord-out solve:', 'DEGENERATE (stage vacuous)' if sol is None else ('unique, 3rd check ok=%s'%sol[2]))
    if sol: print('     vab = (%d, %d)'%(sol[0],sol[1]))
    print('  MUX OUTPUT WIRES: x35346=%d'%(v[35346]%p))
    print('                    x10824=%d'%(v[10824]%p))
    A=(C[10424],C[27436]); B=(C[20930],C[30632])
    cands={"A=(C10424,C27436)":(A[0]%p,A[1]%p),"B=(C20930,C30632)":(B[0]%p,B[1]%p),'A+B (SUM)':((A[0]+B[0])%p,(A[1]+B[1])%p)}
    for orient,(P,Q) in {'coord1=x':(A,B),'coord2=x':((A[1],A[0]),(B[1],B[0]))}.items():
        cx,cy=chordK(P,Q)
        cands['chordK %s'%orient]=(cx,cy) if orient=='coord1=x' else (cy,cx)
    got=(v[35346]%p,v[10824]%p)
    for k,val in cands.items():
        print('     match %-24s : %s'%(k,'YES' if val==got else 'no'))
    # local atom check
    loc=[a for a in E.res if any(w in a for w in ('x36193','x35256','x10424','x27436','x20930','x30632'))]
    bad=[a for a in loc if r[E.residx[a]]%p]
    print('  local atoms: %d, nonzero mod p: %d %s'%(len(loc),len(bad),[a[:60] for a in bad]))

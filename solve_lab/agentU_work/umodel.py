"""Agent-U-successor: slot-pricing model.
Topology from L's full_model.pkl + calib2.pkl (read-only), values over Z (lifted 296-bit),
scored with checker.py, propagated with M's engine (mirror/).
"""
import sys, pickle
sys.set_int_max_str_digits(2000000)
LDIR='/home/user/integer_solver/solve_lab/agentU_work/mirror/L/'
M=pickle.load(open(LDIR+'full_model.pkl','rb'))
C=pickle.load(open(LDIR+'calib2.pkl','rb'))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; live=set(M['live']); link=M['link']
sub=M['sub']; order=M['order']; PIN=M['PIN']; ROOT=M['ROOT']; leafnode=M['leafnode']
ORIENT=C['ORIENT']; perm=C['perm']
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891

XY=pickle.load(open('/home/user/integer_solver/solve_lab/agentU_work/w_xy.pkl','rb'))
LV=pickle.load(open('/home/user/integer_solver/solve_lab/agentU_work/v_leaves.pkl','rb'))
sel2exp=LV['sel2exp']; exp2sel={e:s for s,e in sel2exp.items()}
WZ=pickle.load(open('/home/user/integer_solver/solve_lab/agentU_work/w_z.pkl','rb'))

# lifted per-wire constant + modulus + zvar, keyed by wire id
WCON={}
for e,(w,m,z,Cc) in XY['X'].items(): WCON[w]=(m,z,Cc,'X',e)
for e,(w,m,z,Cc) in XY['Y'].items(): WCON[w]=(m,z,Cc,'Y',e)
# lifted leaf point in PIN wire order, and canonical {'X':..,'Y':..}
LIFT={s:tuple(WCON[w][2] for w in PIN[s][0]) for s in PIN}
AXIS={s:tuple(WCON[w][3] for w in PIN[s][0]) for s in PIN}   # e.g. ('Y','X') or ('X','Y')
LIFTC={s:{WCON[w][3]:WCON[w][2] for w in PIN[s][0]} for s in PIN}
def to_pin(s, canon):
    """canon = {'X':vx,'Y':vy} -> tuple in PIN[s][0] wire order"""
    return tuple(canon[a] for a in AXIS[s])

# parent map + depth
par={}
for n in tree:
    if tree[n]:
        for ch in tree[n]: par[ch]=n
depth={}
def _dep(n):
    if n in depth: return depth[n]
    if tree.get(n) is None: depth[n]=0; return 0
    a,b=tree[n]; d=1+max(_dep(a),_dep(b)); depth[n]=d; return d
for n in order: _dep(n)

LIVELEAF={n:[x for x in sub[n] if x in live] for n in sub}
# the 255 real merge slots
SLOTS=[n for n in order if LIVELEAF[tree[n][0]] and LIVELEAF[tree[n][1]]]

def chord(A,B,o):
    ax,ay,bx,by=A[o],A[1-o],B[o],B[1-o]
    if (bx-ax)%p==0: return None
    l=(by-ay)*pow(bx-ax,p-2,p)%p
    ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    return (ox,oy) if o==0 else (oy,ox)

def buildvals(S, routeval, beta=None, betaval=None):
    """routeval: sel -> 2-tuple (in PIN[sel][0] wire order) used for ROUTING."""
    isl={}; valn={}
    for L in tree:
        if tree[L] is None:
            isl[L]=L in S
            if L in S:
                c=routeval.get(L)
                valn[L]=to_pin(L,c) if c is not None else LIFT.get(L)
            else: valn[L]=None
    for n in order:
        a,b=tree[n]; la,lb=isl[a],isl[b]; isl[n]=la or lb
        def proj(ch,side):
            pm=perm[(n,side)]
            if valn[ch] is None or pm[0] is None or pm[1] is None: return None
            return (valn[ch][pm[0]],valn[ch][pm[1]])
        if n==beta:
            valn[n]=betaval
        elif la and lb:
            o=ORIENT.get(n); PA=proj(a,'va'); PB=proj(b,'vb')
            valn[n]=None if (o is None or o=='DEAD' or PA is None or PB is None) else chord(PA,PB,o)
        elif la: valn[n]=proj(a,'va')
        elif lb: valn[n]=proj(b,'vb')
        else: valn[n]=None
    return isl,valn

def assignment(S, routeval=None, pinval=None, beta=None, betaval=None):
    """Returns dict var->value.  pinval: sel -> 2-tuple written onto the PIN wires."""
    routeval=routeval or {}; pinval=pinval or {}
    isl,valn=buildvals(S,routeval,beta,betaval)
    v={}
    for L in S:
        ws=PIN[L][0]; v[L]=1
        pc=pinval.get(L)
        vals=to_pin(L,pc) if pc is not None else LIFT[L]
        for w,c in zip(ws,vals):
            v[w]=c
            m,z,Cc,ax,e=WCON[w]
            d=c-Cc
            if d: 
                assert d%m==0, 'pin not satisfiable: wire %d'%w
                v[z]=d//m
    for n in order:
        a,b=tree[n]
        for i,d in enumerate(OUT[n]):
            for side,ch in (('va',a),('vb',b)):
                if tree[ch] is not None:
                    pmi=perm[(n,side)][i]
                    v[d[side]]=(valn[ch][pmi] if (isl[ch] and pmi is not None and valn[ch] is not None) else 0)
            v[d['vab']]=(valn[n][i] if (valn[n] is not None and (n==beta or (isl[a] and isl[b]))) else 0)
    return v,isl,valn

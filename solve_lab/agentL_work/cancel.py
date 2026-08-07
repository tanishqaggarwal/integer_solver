"""Cancellation search: generalise the deliverable's actual cut and price it EXACTLY.

Deliverable's real mechanism (measured, delivsite.py/rootcheck.py):
  ON = {24601, 2081}.  The 2081 side is corrupted so that at the ROOT both inputs are the
  SAME value (leaf 24601's).  With the two inputs equal the root's own 3 stage checks are
  satisfiable for an arbitrary output, so the root's vab wires are set straight to TARGET and
  everything above closes.  Cost = forging that equality: 2 guards at x27994 + 2 slot links
  at x4971.va = 4 atoms = 7 equations.
"""
import sys, os, json, pickle, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import checker as CK
src=open('/home/user/integer_solver/solve_lab/agentL_work/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
TGT=tuple(pickle.load(open('target.pkl','rb')))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; sub=M['sub']; ROOT=M['ROOT']; liveset=set(M['live'])
parent={}; side_of={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n; side_of[ch]=s
def chordf(A,B,o):
    ax,ay,bx,by=A[o],A[1-o],B[o],B[1-o]
    d=(bx-ax)%p
    if d==0: return None
    l=(by-ay)*pow(d,p-2,p)%p
    ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    return (ox,oy) if o==0 else (oy,ox)
def invchord(O,A,o):
    ax,ay,ox,oy=A[o],A[1-o],O[o],O[1-o]
    d=(ax-ox)%p
    if d==0: return None
    l=(oy+ay)*pow(d,p-2,p)%p
    bx=(l*l-ax-ox-K)%p; by=(ay+l*(bx-ax))%p
    return (bx,by) if o==0 else (by,bx)
print('loading checker equations...',flush=True)
t0=time.time(); CODES,_=CK.load_equations(); print('  %d eqs in %.0fs'%(len(CODES),time.time()-t0),flush=True)
NS_V=[0]*CK.NVARS
def exact_fail(vv):
    v=[0]*CK.NVARS
    n=min(len(vv),CK.NVARS)
    v[:n]=vv[:n]
    ns={'v':v,'__builtins__':{}}
    return sum(1 for c in CODES if eval(c,ns)!=0)
if __name__=='__main__':
    D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
    vd=[0]*CK.NVARS
    for k,x in D.items(): vd[int(k[2:])]=int(x)
    print('sanity: exact_fail(deliverable) =',exact_fail(vd),'(expect 7)',flush=True)
    va=json.load(open('assign_L1.json'))
    v2=[0]*CK.NVARS
    for k,x in va.items(): v2[int(k[2:])]=int(x)
    print('sanity: exact_fail(my assign_L1) =',exact_fail(v2),'(expect 15)',flush=True)

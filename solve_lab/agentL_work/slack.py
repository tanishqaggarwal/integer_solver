"""Is anything forcing x_4116 and its sibling shared factors to zero?"""
import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentF_work')
from fwd import Engine,NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
Hd=pickle.load(open('handles.pkl','rb'))
Z=set(Hd.get('Z',[])) if 'Z' in Hd else set()
handle=set(Hd['handle']); value=set(Hd['value'])
SHARED=[4116,16153,1962,12682,19049,15616]
OTHERF=[22163,10858,14393]
SLACK=[36780,11630,34243]
def show(v,tag):
    d=defrhs.get(v)
    print('  x%-6d %-8s free=%-5s  def=%s'%(v,tag,v not in defrhs, node_str(d)[:90] if d else '-'))
    if d is None: return
print('--- the six shared factors Q flagged ---')
for v in SHARED: show(v,'SHARED')
print('--- the other factor in each slack product ---')
for v in OTHERF: show(v,'FACTOR')
print('--- the slack wires themselves ---')
for v in SLACK: show(v,'SLACK')
print()
print('p =',p)
# evaluate everything on the all-zero assignment: constants show up immediately
vv=[0]*NV
E.run(vv)
V=E.vals(vv) if hasattr(E,'vals') else None
# fall back: evaluate defs directly
def ev(v,seen=None):
    if v not in defrhs: return None
    r=defrhs[v]
    def e(n):
        if n[0]=='c': return n[1]
        if n[0]=='v': return ev(n[1])
        a=e(n[1]); b=e(n[2])
        if a is None or b is None: return None
        return a+b if n[0]=='+' else (a-b if n[0]=='-' else a*b)
    return e(r) if r[0] in '+-*' else (r[1] if r[0]=='c' else ev(r[1]))
sys.setrecursionlimit(100000)
print('--- CONSTANT VALUE of each (evaluated from definitions alone, no free vars) ---')
for v in SHARED+OTHERF+SLACK:
    val=ev(v)
    print('  x%-6d = %s   == p ? %s   ==0 ? %s'%(v,str(val)[:80] if val is not None else 'depends on free vars',
          val==p if val is not None else '-', val==0 if val is not None else '-'))

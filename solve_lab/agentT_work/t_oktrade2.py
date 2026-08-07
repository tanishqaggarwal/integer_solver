#!/usr/bin/env python3
"""AUDIT T28b -- the crux.  T28 found 7 knobs that move a failing row while leaving S fixed, so
O's 'every purchase costs exactly eq8680' is NOT forced by 'all knobs move S'.  But MOVING a row
is not BUYING it.  Test whether the S-preserving subspace can actually zero a failing row."""
import os,sys,json,pickle,re,collections
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentI_work'))
import checker as CK
from intsolve import solve_int
NV=38748
B=json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json')))
v0=[0]*NV
for k,val in B.items(): v0[int(k[2:])]=int(val)
lhs=open('/home/user/integer_solver/EQUATIONS.txt').read().split('\n')[8680].rsplit('=',1)[0].strip()
def so(s):
    while s.startswith('(') and s.endswith(')'):
        dd=0; ok=True
        for k,ch in enumerate(s):
            if ch=='(': dd+=1
            elif ch==')':
                dd-=1
                if dd==0 and k!=len(s)-1: ok=False; break
        if not ok: break
        s=s[1:-1]
    return s
def sp(s,op):
    s=so(s); out=[]; dd=0; cur=''
    for ch in s:
        if ch=='(': dd+=1
        elif ch==')': dd-=1
        if dd==0 and ch==op: out.append(cur); cur=''; continue
        cur+=ch
    out.append(cur); return out
S=so(sp(so(sp(lhs,'*')[0]),'*')[0])
code_S=compile(re.sub(r'x_(\d+)',r'v[\1]',S),'<S>','eval')
codes,_=CK.load_equations()
FAILS=[12231,12270,12350,14584,18673,22044,29125]
SPRES=[1329,7068,8731,9118,9413,10903,17325]     # move a failing row, leave S fixed (T28)
print('S-preserving knobs: %s'%SPRES)
base_e={e:eval(codes[e],{'v':v0,'__builtins__':{}}) for e in FAILS}
baseS=eval(code_S,{'v':v0,'__builtins__':{}})
# affinity + jacobian over the S-preserving knobs
J={e:[] for e in FAILS}; JS=[]
aff=True
for u in SPRES:
    d1=[];
    for t in (1,2):
        vv=list(v0); vv[u]=v0[u]+t
        ns={'v':vv,'__builtins__':{}}
        d1.append(({e:eval(codes[e],ns)-base_e[e] for e in FAILS}, eval(code_S,ns)-baseS))
    for e in FAILS:
        if d1[1][0][e]!=2*d1[0][0][e]: aff=False
    if d1[1][1]!=2*d1[0][1]: aff=False
    for e in FAILS: J[e].append(d1[0][0][e])
    JS.append(d1[0][1])
print('all rows + S affine in these knobs? %s'%aff)
print('dS/dknob over the S-preserving set: %s'%JS)
print('\ncan the S-preserving subspace ZERO a failing row (integer solve)?')
for e in FAILS:
    x=solve_int([J[e]],[-base_e[e]])
    # also require S stay 0
    x2=solve_int([J[e],JS],[-base_e[e],0])
    print('   eq%-6d  row alone: %-9s   row AND S=0 kept: %s'%(
        e,'SOLVABLE' if x is not None else 'no',
        'SOLVABLE' if x2 is not None else 'no'))
print('\ninterpretation:')
print('  if no row is solvable with S held at 0, O\'s uniformity is explained after all --')
print('  not because every knob moves S, but because the S-preserving directions cannot')
print('  reach any row\'s target.  That is a genuine search result, not a restatement.')

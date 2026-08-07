#!/usr/bin/env python3
"""AUDIT T28 -- agent O's seven-way trade.  Two angles, both in F's certified-faithful parse.

ANGLE 1: is K what it says it is?  O: K = (15 free inputs reaching a nonzero region atom)
         u (26 carriers of S) = 34.  Note 15+26 = 41, so O implies an overlap of exactly 7.
ANGLE 2: is 'every purchase costs exactly eq8680' ONE fact seen seven times, or seven
         independent measurements that agree?  eq8680 = S^4 (my 8th pass) and S = 0 at the
         witness, so a purchase costs eq8680 iff it moves S off zero.  If EVERY knob that can
         move a failing equation also moves S, the uniformity is structurally forced -- one fact.
         If some knob moves a failing row while leaving S fixed, the seven are independent."""
import os,sys,json,pickle,collections,re
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from fwd import compile_node
from circ2 import vars_of
import checker as CK, engine as E
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
NV=38748
B=json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json')))
v0=[0]*NV
for k,val in B.items(): v0[int(k[2:])]=int(val)
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<at>','exec')
r=[0]*len(names); exec(prog,{'v':v0,'r':r,'__builtins__':{}})
REGION=[names[i] for i in range(len(names)) if r[i]]
print('region atoms (nonzero at the witness): %d'%len(REGION))
# --- S, from the raw file (8th pass) ---
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
Sv=eval(code_S,{'v':v0,'__builtins__':{}})
print('S at the witness = %s  (eq8680 satisfied: %s)'%(Sv,Sv==0))
Scar=sorted({int(m) for m in re.findall(r'x_(\d+)',S)})
Sfree=[u for u in Scar if E.definer[u] is None]
print('\nANGLE 1')
print('   S carriers: %d variables, of which FREE: %d   (O says 26)'%(len(Scar),len(Sfree)))
# free inputs reaching a region atom, through the definition DAG
alldef={}
for a in names:
    m=re.match(r'^\(x(\d+)-(.+)\)$',a.replace(' ',''))
    if m:
        vv=int(m.group(1))
        if vv not in alldef: alldef[vv]=vars_of(atoms[a])-{vv}
def freein(start):
    seen=set(); st=list(start); out=set()
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        if u in alldef:
            for w in alldef[u]: st.append(w)
        else: out.add(u)
    return out
reg=set()
for a in REGION: reg|=vars_of(atoms[a])
Rfree=sorted(freein(reg))
print('   free inputs reaching a region atom: %d   (O says 15)'%len(Rfree))
K=sorted(set(Rfree)|set(Sfree))
print('   |K| = |union| = %d   (O says 34)   overlap = %d'%(len(K),len(set(Rfree)&set(Sfree))))
print('\nANGLE 2  -- does every knob that moves a failing row also move S?')
codes,_=CK.load_equations()
FAILS=[12231,12270,12350,14584,18673,22044,29125]
eqcodes={e:codes[e] for e in FAILS}
base={e:eval(codes[e],{'v':v0,'__builtins__':{}}) for e in FAILS}
print('   %-9s %-9s %s'%('knob','moves S?','failing rows it moves'))
movesS=0; movesRowNotS=[]
for u in K:
    vv=list(v0); vv[u]=v0[u]+1
    ns={'v':vv,'__builtins__':{}}
    dS = eval(code_S,ns)-Sv
    mv=[e for e in FAILS if eval(eqcodes[e],ns)!=base[e]]
    if dS: movesS+=1
    if mv and not dS: movesRowNotS.append((u,mv))
    if mv or dS:
        print('   x%-8d %-9s %s'%(u,'YES' if dS else 'no',mv if mv else '-'))
print('\n   knobs in K that move S: %d of %d'%(movesS,len(K)))
print('   knobs that move a FAILING row but leave S fixed: %d  %s'%(len(movesRowNotS),movesRowNotS[:5]))
if not movesRowNotS:
    print('   -> EVERY knob in K that can touch a failing row also moves S.  Since S = 0 is ONE')
    print('      linear constraint and eq8680 is its only equation, "every purchase costs exactly')
    print('      eq8680" is STRUCTURALLY FORCED: one fact seen seven times, not seven agreeing')
    print('      measurements.  The 1-for-1 trade cannot be leveraged for the same reason.')
else:
    print('   -> some knob moves a failing row without moving S: the seven are NOT one fact.')

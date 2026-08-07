#!/usr/bin/env python3
"""Q-9a: re-derive the leaf pins DIRECTLY from EQUATIONS.txt (no dependence on any other agent).
Pin atoms have the shape  (x_g)*((x_w)-(BIGCONST)) : when selector x_g = 1 the coordinate wire
x_w is forced to BIGCONST."""
import re,json,collections
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
cs = K*pow(3,p-2,p)%p
PIN=re.compile(r'\(x_(\d+)\)\)?\*\(\(x_(\d+)\)-\((\d{40,})\)\)')
txt=open('/home/user/integer_solver/EQUATIONS.txt').read()
pins=collections.defaultdict(dict)
for g,w,c in PIN.findall(txt):
    pins[int(g)][int(w)]=int(c)
print('selectors carrying big pins:',len(pins))
cnt=collections.Counter(len(v) for v in pins.values())
print('pins per selector:',dict(cnt))
json.dump({str(g):{str(w):str(c) for w,c in v.items()} for g,v in pins.items()},open('qpins.json','w'))
# each selector with 2 pins => a leaf point.  Which of the two wires is X and which is Y?
a,b = int(json.load(open('curve.json'))['a']), int(json.load(open('curve.json'))['b'])
def oncur(X,Y): return (Y*Y-pow(X,3,p)-a*X-b)%p==0
good=0; bad=0; leaf={}
for g,v in pins.items():
    if len(v)!=2: continue
    (w1,c1),(w2,c2)=sorted(v.items())
    X1,Y1=(c1+cs)%p, c2%p
    X2,Y2=(c2+cs)%p, c1%p
    if oncur(X1,Y1): leaf[g]=((X1,Y1),(w1,w2)); good+=1
    elif oncur(X2,Y2): leaf[g]=((X2,Y2),(w2,w1)); good+=1
    else: bad+=1
print('two-pin selectors giving an on-curve point: %d   off-curve: %d'%(good,bad))
json.dump({str(g):[str(v[0][0]),str(v[0][1]),v[1][0],v[1][1]] for g,v in leaf.items()},open('qleaf.json','w'))

"""Find the load-pin atoms  bit*(x_B - HUGE) - s*x_C  and group them by selector bit.
These are the structurally-gated leaf points."""
import sys, collections, json; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
pins=[]
for a in range(L.NA):
    Pp=L.polys[a]
    big=[(m,c) for m,c in Pp.items() if abs(c)>=10**40]
    if not big: continue
    # want: exactly one huge coefficient, on a degree-1 monomial (the gate bit) OR deg1
    # general shape: sum of terms; identify the monomial carrying the huge constant
    ms=sorted(Pp.items(), key=lambda kv:-abs(kv[1]))
    pins.append((a,Pp))
print('atoms containing a >=40-digit literal: %d'%len(pins))
shapes=collections.Counter()
for a,Pp in pins:
    sig=tuple(sorted((len(m), 'H' if abs(c)>=10**40 else ('p' if abs(c)==P else 's')) for m,c in Pp.items()))
    shapes[sig]+=1
for s,c in shapes.most_common(12): print('   %-5d %s'%(c,s))
# focus: shape with exactly one huge term of degree 1 (the constant times a bit)
sel=[]
for a,Pp in pins:
    hs=[(m,c) for m,c in Pp.items() if abs(c)>=10**40]
    if len(hs)!=1: continue
    m,c=hs[0]
    if len(m)!=1: continue
    bit=m[0]
    rest=[(mm,cc) for mm,cc in Pp.items() if (mm,cc)!=(m,c)]
    sel.append((a,bit,c,rest))
print('\npins with exactly ONE huge term on a single variable: %d'%len(sel))
bybit=collections.defaultdict(list)
for a,bit,c,rest in sel: bybit[bit].append((a,c,rest))
mult=collections.Counter(len(x) for x in bybit.values())
print('pins per gating variable histogram:',sorted(mult.items()))
print('gating variables with EXACTLY 2 pins: %d'%sum(1 for x in bybit.values() if len(x)==2))
print('distinct gating variables: %d'%len(bybit))
ex=[b for b,x in bybit.items() if len(x)==2][:4]
for b in ex:
    print('  gate x%d (value %s):'%(b,v[b]))
    for a,c,rest in bybit[b]:
        rs=' '.join(('%+d'%cc if abs(cc)<10**12 else '%+d~%dd'%(cc%1000,len(str(abs(cc)))))+''.join('*x%d'%u for u in mm) for mm,cc in sorted(rest))
        print('     a%-6d  HUGE=%d  rest= %s'%(a,c,rs))
json.dump({str(b):[[a,str(c)] for a,c,_ in x] for b,x in bybit.items()},
          open('/home/user/integer_solver/solve_lab/agentA_work/pins.json','w'))

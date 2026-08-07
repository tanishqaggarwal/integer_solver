"""Measure the addition-law offset K directly from the instance, then test the
'pin constants on y^2 = (x+s)^3 + b' hypothesis with that K (and by solving for s)."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
X1,Y1,X2,Y2,X3,Y3=12186,16742,14853,24908,22162,30213
for t in [35389,6671]:
    a=L.definer.get(t)
    print('x%d definer atom = %s'%(t,a))
    if a is not None:
        Pp=L.polys[a]
        print('   %d terms, vars=%s'%(len(Pp),sorted(L.avars[a])[:25]))
# measure A as a function of x3 by finite differences over Z: A = c2*(x3) + c1  when x1,x2,y1,y2 fixed
def val(t,w): return w[t]
def probe(target, knob, deltas):
    out=[]
    for d in deltas:
        w=list(v); w[knob]=w[knob]+d
        out.append(L.evalpoly(L.polys[L.definer[target]],w) if target in L.definer else None)
    return out
# Instead: measure the identity directly with the CURRENT values
x1,y1,x2,y2,x3,y3=[v[t]%P for t in (X1,Y1,X2,Y2,X3,Y3)]
print('\nx1==x2 ? %s   y1==y2 ? %s'%(x1==x2,y1==y2))
# K from A=0: (x2-x1)^2*(x3+x1+x2+K) = (y2-y1)^2  -> undefined when x1==x2
# use AG_39013 where x1 != x2
w=L.load('/home/user/integer_solver/solve_lab/s10/AG_39013.json')
x1,y1,x2,y2,x3,y3=[w[t]%P for t in (X1,Y1,X2,Y2,X3,Y3)]
lam=(y2-y1)*pow(x2-x1,-1,P)%P
K=(lam*lam-x1-x2-x3)%P
print('AG_39013: implied K = lam^2 - x1 - x2 - x3 = %d'%K)
print('   K/3 mod p = %d'%(K*pow(3,-1,P)%P))
# test: do the literals pair on y^2=(x+s)^3+b for s = K/3 ?
lits=set()
for a in range(L.NA):
    for m_,c in L.polys[a].items():
        if abs(c)>=10**40: lits.add(abs(c)%P)
lits=sorted(lits)
sq=collections.defaultdict(list)
for Y in lits: sq[pow(Y,2,P)].append(Y)
b1=(v[Y1]*v[Y1]-pow(v[X1],3,P))%P
for name,s in [('0',0),('K',K),('K/3',K*pow(3,-1,P)%P),('-K/3',(-K)*pow(3,-1,P)%P)]:
    for bname,b in [('b1',b1),('7',7)]:
        c=0
        for X in lits:
            if ((pow((X+s)%P,3,P)+b)%P) in sq: c+=1
        if c: print('   s=%s b=%s -> %d literal pairs'%(name,bname,c))
print('done')

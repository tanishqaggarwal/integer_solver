"""INDEPENDENT test of the 'the pin constants lie on y^2 = x^3 + b' hypothesis.
Collect every large literal coefficient in the atom polynomials, reduce mod p, and look
for a b such that many ordered pairs (X,Y) satisfy Y^2 = X^3 + b (mod p)."""
import sys, collections, json; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
lits=collections.Counter()
for a in range(L.NA):
    for m,c in L.polys[a].items():
        if abs(c)>=10**40: lits[abs(c)]+=1
print('large literals (>=40 digits): %d distinct, %d occurrences'%(len(lits),sum(lits.values())))
bits=collections.Counter(len(bin(x))-2 for x in lits)
print('bit-length histogram:',sorted(bits.items()))
Ls=sorted(lits)
R=[x%P for x in Ls]
print('distinct residues mod p: %d'%len(set(R)))
# b = Y^2 - X^3 for all ordered pairs
bc=collections.Counter()
for X in set(R):
    x3=pow(X,3,P)
    for Y in set(R):
        bc[(pow(Y,2,P)-x3)%P]+=1
top=bc.most_common(8)
print('most common b = Y^2 - X^3 over all ordered pairs of literal residues:')
for b,cnt in top:
    print('   count=%-6d b=%d%s'%(cnt,b,'   <== b == 7 (secp256k1)' if b==7 else ''))
n=len(set(R))
print('expected count for a random b ~ %.2f  (n^2/p is ~0)'%(n*n/ max(len(bc),1)))
json.dump({'n_literals':len(lits),'top_b':[[str(b),c] for b,c in top]},
          open('/home/user/integer_solver/solve_lab/agentA_work/curve_test.json','w'))

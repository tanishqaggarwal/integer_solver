"""Close the two mod-p congruences at mod9118_0 by shifting the free inputs x14853, x14623."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s10')
import tools as T, ad
P=env.P
ORDER=ad.ORDER
def fwd(v,rounds=6):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(L.definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
def score(v):
    av=L.all_atom_values(v); fe=L.failing_eqs(av)
    return L.NEQ-len(fe), fe, [a for a in range(L.NA) if av[a]]
path=sys.argv[1]
v=L.load(path)
s0,fe0,nz0=score(v); print('load %s -> %d, nz=%s'%(path.split('/')[-1],s0,nz0))
w=list(v); fwd(w); s1,fe1,nz1=score(w); print('after canonical fwd -> %d, nz=%s'%(s1,nz1))
base = w if s1>=s0 else list(v)
print('using base score',score(base)[0])
# congruence A: a21617 = c1*x14623 + c2*x27522 + c3*x36864  (mod p)
# congruence B: a29539 = d1*x1308 + d2*x14853 + d3*x29967
def lincoef(a,u,v):
    c=0
    for m,cc in L.polys[a].items():
        if u in m:
            t=cc
            for w2 in m:
                if w2!=u: t*=v[w2]
            c+=t
    return c
for (atom,u) in [(21617,14623),(29539,14853)]:
    av=L.all_atom_values(base)
    c=lincoef(atom,u,base)%P
    val=av[atom]%P
    if c==0: print('  a%d: coeff of x%d is 0 mod p -- cannot fix'%(atom,u)); continue
    d=(-val)*pow(c,-1,P)%P
    print('  a%d val mod p=%d, coeff on x%d=%d, shift delta=%d'%(atom,val,u,c,d))
    base[u]=base[u]+d
    av2=L.all_atom_values(base)
    print('     -> a%d mod p now %d'%(atom,av2[atom]%P))
s2,fe2,nz2=score(base); print('after shifts (no fwd) -> %d nz=%s'%(s2,nz2))
fwd(base)
s3,fe3,nz3=score(base); print('after fwd -> %d nz=%s failing=%s'%(s3,fe3 and len(fe3),nz3))
print('SCORE %d'%s3)
if s3>=39026:
    out='/home/user/integer_solver/solve_lab/agentA_work/A_exp1_%d.json'%s3
    json.dump({str(i):str(base[i]) for i in range(L.NVARS)},open(out,'w')); print('saved',out)

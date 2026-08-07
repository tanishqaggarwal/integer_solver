"""Is there a chain of EC points on ONE curve inside the state?
Collect all distinct large variable values mod p, and count how many ordered pairs (X,Y)
satisfy Y^2 = X^3 + b for each b.  A doubling ladder over 256 steps would show a b with
>= ~256 pairs.  Also report the maximum multiplicity of any b."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
v0=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
b1=(v0[16742]*v0[16742]-pow(v0[12186],3,P))%P
for path in sys.argv[1:]:
    v=L.load(path)
    vals=sorted({x%P for x in v if abs(x)>2**200})
    print('%s : %d distinct large variable residues'%(path.split('/')[-1],len(vals)))
    sq=collections.defaultdict(list)
    for Y in vals: sq[pow(Y,2,P)].append(Y)
    for name,b in [('b1',b1),('7',7)]:
        hits=[(X,sq[(pow(X,3,P)+b)%P][0]) for X in vals if ((pow(X,3,P)+b)%P) in sq]
        print('   pairs on y^2=x^3+%s : %d %s'%(name,len(hits),hits[:2]))
    # global: max multiplicity of any b over all ordered pairs (only if small enough)
    if len(vals)<=4000:
        bc=collections.Counter()
        cubes=[pow(X,3,P) for X in vals]
        squares=[pow(Y,2,P) for Y in vals]
        for c in cubes:
            for s in squares: bc[(s-c)%P]+=1
        mc=bc.most_common(4)
        print('   max b-multiplicity over %d^2 pairs: %s'%(len(vals),[(str(b)[:20]+'..',c) for b,c in mc]))

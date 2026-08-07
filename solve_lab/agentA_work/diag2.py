import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
v=L.load(sys.argv[1]); av=L.all_atom_values(v)
for a in [int(x) for x in sys.argv[2:]]:
    print('=== a%d  val mod p = %d ==='%(a, av[a]%P))
    for m,c in sorted(L.polys[a].items()):
        t=c
        for u in m: t*=v[u]
        info=[]
        for u in m:
            n=len(L.var_atoms[u]); d='free' if u not in L.definer else 'G(a%d)'%L.definer[u]
            info.append('x%d[%s,%datoms,%s]'%(u,d,n,'p' if v[u]==P else ('0' if v[u]==0 else ('%d'%v[u] if abs(v[u])<10**6 else '%dd'%len(str(abs(v[u])))))))
        print('   %s   termmodp=%d' % (' * '.join(info), t%P))

"""What does the instance ASSERT?  Look at the 256 selector bits' values and hunt for the
target point Q of the scalar multiplication among the pinned/literal points."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
c=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/curve_final.json'))
sh=int(c['shift']); a2=int(c['a2']); a4=int(c['a4']); a6=int(c['a6']); B=int(c['B'])
d=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/weier.json'))
NSEC=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
def onc(x,y): return (pow(y,2,P)-pow(x,3,P)-a2*pow(x,2,P)-a4*x-a6)%P==0
for path in ['/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json',
             '/home/user/integer_solver/solve_lab/s10/AG_39013.json']:
    v=L.load(path)
    bits=sorted(int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentA_work/pins.json')))
    on=[b for b in bits if v[b]!=0]
    print('%s : selector bits set = %d of %d  %s'%(path.split('/')[-1],len(on),len(bits),on[:12]))
    # any pair of VARIABLES forming a point on the curve?
    vals={u:v[u]%P for u in range(L.NVARS) if abs(v[u])>2**200}
    pts=[]
    items=list(vals.items())
    sq={}
    for u,y in items: sq.setdefault(y,[]).append(u)
    for u,x in items:
        rhs=(pow(x,3,P)+a2*pow(x,2,P)+a4*x+a6)%P
        for w,y in items:
            if (y*y-rhs)%P==0: pts.append((u,w,x,y))
    print('   variable pairs (u,w) with (v[u],v[w]) ON the curve: %d'%len(pts))
    for u,w,x,y in pts[:8]:
        print('      x%-6d , x%-6d   x=%d'%(u,w,x))
    # the pinned point
    x1,y1=v[12186]%P,v[16742]%P
    print('   (x12186,x16742) on curve: %s'%onc(x1,y1))
    # literal points not equal to a 2^i*G leaf?  check the 296-bit literal residues

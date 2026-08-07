import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close import closure, evalat, CHK
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
y1=75252414949035169209861377566271802862155171514703932453490974720978433969307
x2=115654104920822720384433381653801116189600145117015552377208942221379503483440
for tag,extra in [('stall',{}),('AB0',{16742:y1,14853:x2})]:
    s=dict(BASE); s.update(extra)
    seeds,v,sc,nz=closure(s,verbose=False)
    print('=== %s score=%d nz=%d'%(tag,sc,len(nz)))
    print('   x_11150=%d x_25739=%d x_37758=%d'%(v[11150]%P,v[25739]%P,v[37758]%P))
    for a in nz:
        print('   a%-6d eqs=%-3d val=%s'%(a,len(L.atom2eq.get(a,{})),str(evalat(a,v))[:36]))
        print('       ',L.atom_src[a][:230])

import sys, json, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
def run(extra):
    c=dict(BASE); c.update(extra)
    sc,v,nz=closure2(c)
    return sc,v,nz
random.seed(3)
base_sc,bv,bnz=run({})
def rep(v,tag,sc,nz):
    print('[%s] sc=%d nz=%d  y1=%d.. x2=%d.. A=%d.. B=%d..'%(tag,sc,len(nz),
        v[16742]%P%10**12,v[14853]%P%10**12,v[35389]%P%10**12,v[6671]%P%10**12))
rep(bv,'base',base_sc,bnz)
for name,u in [('x8778',8778),('x6418',6418)]:
    for t in [1, 12345678901234567890, random.randrange(P)]:
        sc,v,nz=run({u:t})
        rep(v,'%s=%s'%(name,str(t)[:12]),sc,nz)
        print('     x_19083=%d.. x_1308=%d.. 11150=%d.. 25739=%d.. 37758=%d..'%(v[19083]%P%10**10,v[1308]%P%10**10,v[11150]%P%10**10,v[25739]%P%10**10,v[37758]%P%10**10))

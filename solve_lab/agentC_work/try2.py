import sys, json, re
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
CHK=[a for a in range(L.NA) if a not in L.atom_out]
def run(seeds,tag,show=0):
    v=[0]*L.NVARS
    for k,val in seeds.items(): v[k]=val
    forward(v)
    ov=[k for k,val in seeds.items() if v[k]!=val]
    av=L.all_atom_values(v)
    f=L.failing_eqs(av)
    nz=[a for a in CHK if av[a]!=0]
    print('[%s] score=%d nz=%d overwritten=%s'%(tag,L.NEQ-len(f),len(nz),ov))
    for a in nz[:show]:
        print('   a%-6d val=%s'%(a,str(av[a])[:50])); print('      ',L.atom_src[a][:220])
    return v,av,nz,f
S={542:1, 91:1, 22162:K2, 30213:K1,
   13153:26066848117641554604733072514062377518640010461224574295161576128745599678073481099486949,
   20386:123605134436201557989718624451460897015798857659190965912156591009360284284149351426198201,
   35344:22742974578460079609867846888465534794042730731697871517291778569302981207665435077363969,
   23210:12305317472058964629688036168955847079074597599772814644736410753792807337224143169444663}
v,av,nz,f=run(S,'pins-set',show=20)
json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0}, open('/home/user/integer_solver/solve_lab/agentC_work/try2.json','w'))

import pickle, sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M=pickle.load(open('full_model.pkl','rb')); C2=pickle.load(open('calib2.pkl','rb'))
link=M['link']; OUT=M['OUT']; ROOT=M['ROOT']; tree=M['tree']; live=M['live']; sub=M['sub']
T1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
T2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
outs=[d['out'] for d in OUT[ROOT]]
print('root out wires',outs,' link[24468]=',link.get(24468),' link[18956]=',link.get(18956))
tgt=[None,None]
for w,T in ((24468,T1),(18956,T2)):
    z=link.get(w)
    if z in outs: tgt[outs.index(z)]=T%p
print('TARGET (in root coord order):',tgt)
a,b=tree[ROOT]
print('root children x%d x%d ; live leaves %d | %d'%(a,b,len([x for x in sub[a] if x in set(live)]),len([x for x in sub[b] if x in set(live)])))
# depth profile of live-leaf counts
ls=set(live)
prof=sorted((len([x for x in sub[n] if x in ls]) for n in M['NODE']))
print('nodes with live-leaf support <=24:',sum(1 for x in prof if x<=24),'of',len(prof))
print('largest supports',prof[-12:])
pickle.dump(tgt,open('target.pkl','wb'))

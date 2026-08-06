import heal_harness as H
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
def show(idx):
    print(f"x_{idx} = {H.val[idx]%p}   free={idx in H.freeinp}   ancfree={sorted(H.anc.get(idx,set()))[:8] if idx not in H.freeinp else '-'}")
for idx in [15298,11150,25739,29804,37758,35605,4007, 14853,12186,24908,16742,29322,3558,35389,6671]:
    show(idx)
print()
# the 4 atoms residuals
def atomval(poly): 
    s=0
    for m,c in poly:
        t=c%p
        for v in m: t=(t*H.val[v])%p
        s=(s+t)%p
    return s
A={18081:[((15298,11150),1)],
   18084:[((15298,25739),1),((29804,),-6672769)],
   29377:[((15298,37758),537773),((35605,),-1)],
   35321:[((15298,11150),1),((4007,),1)]}
for i,pol in A.items():
    print(f"atom {i} = {atomval(pol)}")
# check MUX val
print("\nx_15298 value (full int):", H.val[15298])

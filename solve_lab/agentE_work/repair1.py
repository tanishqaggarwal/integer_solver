import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
def show(v,tag,av):
    print(f"{tag}: a=x7715={v[7715]} b=x34554={v[34554]} x13913={str(v[13913])[:30]} x13682={str(v[13682])[:30]} "
          f"x12186={str(v[12186])[:20]} x14853={str(v[14853])[:20]} x24530={str(v[24530])[:25]} x24908={str(v[24908])[:25]} bad={sorted(av)}")
for bit,label in [(22106,'a=1'),(5090,'b=1')]:
    s={18956:C, bit:1}
    v=H.forward(s); ff,av=H.eqfails(v); show(v,label+' step0 fails=%d'%len(ff),av)
    # iterate: pick the mux input to match
    for it in range(6):
        if v[7715]==1 and v[34554]==0: s[12186]=v[13682]
        elif v[34554]==1 and v[7715]==0: s[14853]=v[13682]
        v=H.forward(s); ff,av=H.eqfails(v); show(v,f'{label} it{it} fails={len(ff)}',av)
    json.dump({f"x_{i}":v[i] for i in range(H.NV) if v[i]!=0}, open(f'rep_{label}.json','w'))

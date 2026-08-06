import heal_harness as H
import json, glob
p=H.p
apoly={18081:[((15298,11150),1)],18084:[((15298,25739),1),((29804,),-6672769)],
       29377:[((15298,37758),537773),((35605,),-1)],35321:[((15298,11150),1),((4007,),1)]}
def atomval(items):
    s=0
    for m,c in items:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
# selector-ish bits from cert pin atoms
selectors=[18022,490,34309,24774,1931,3106,11368,22562,28005,16348,29759,11681,24501]
files=sorted(glob.glob('best/*.json'))
print(f"{'file':45s} nfail x_15298 x_11150%p(short) atom18081 selbits")
for fn in files:
    try:
        d=H.loadd(fn)
    except Exception as e:
        print(fn,"load err",e); continue
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    F=H.fails()
    x15=H.val[15298]%p
    x11=H.val[11150]%p
    a18=atomval(apoly[18081])
    sel=''.join(str(H.val[s]%p if H.val[s]%p<10 else '?') for s in selectors)
    print(f"{fn.split('/')[-1]:45s} {len(F):4d}  {x15:1d}   ...{str(x11)[-8:]}  {'0' if a18==0 else 'NZ':2s}  {sel}")

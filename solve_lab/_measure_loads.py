import heal_harness as H, json, time
from collections import defaultdict
p=H.p
pins=json.load(open('pinrec.json'))
bysel=defaultdict(list)
for rec in pins:
    atomidx,sel,tgt,const,coef,handle=rec
    bysel[sel].append((tgt,const%p,handle))
bits=sorted(bysel)
LOADS=[11150,25739,37758]
# baseline all-zero
for v in H.freeinp: H.val[v]=0
H.forward()
base=[H.val[n]%p for n in LOADS]
print("baseline loads:",base, "mux",H.val[15298]%p)
# single-bit response
resp={}   # bit -> (L1,L2,L3) mod p when only this bit active
changed=[]
for b in bits:
    saved={b:H.val[b]}; H.val[b]=1
    for tgt,const,handle in bysel[b]:
        saved[tgt]=H.val[tgt]; H.val[tgt]=const
    H.forward()
    lv=tuple(H.val[n]%p for n in LOADS)
    mv=H.val[15298]%p
    resp[b]=(lv,mv)
    if lv!=tuple(base) or mv!=0: changed.append(b)
    for v,x in saved.items(): H.val[v]=x
H.forward()  # restore
# coefficients c_{k,j} = resp - base
coef={}
for b in bits:
    lv,mv=resp[b]
    coef[b]=tuple((lv[k]-base[k])%p for k in range(3))
nzL=[b for b in bits if any(coef[b][k]!=0 for k in range(3))]
nzMux=[b for b in bits if resp[b][1]!=0]
print("bits changing any load:",len(nzL))
print("bits changing mux=1:",len(nzMux))
# among load-changing bits, breakdown per load
for k,name in enumerate(['L1','L2','L3']):
    s=[b for b in bits if coef[b][k]!=0]
    print(f"  {name}: {len(s)} bits affect it")
json.dump({'base':base,'coef':{str(b):list(coef[b]) for b in bits},
           'mux':{str(b):resp[b][1] for b in bits}},
          open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/loadcoef.json','w'))
print("saved loadcoef.json")
# print a few nonzero coefs
for b in nzL[:6]:
    print("bit",b,"coef",[str(c)[:16] for c in coef[b]])

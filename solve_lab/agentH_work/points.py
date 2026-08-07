"""Extract the 256 leaf points P_b = (H1_b, H2_b) from the bit load pins."""
import ev, pickle, json, os
from fast import chk, csup, inv
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
polys=pickle.load(open(os.path.join(HERE,'polys.pkl'),'rb'))
BITS=json.load(open(os.path.join(HERE,'bits.json')))
allbits=sorted(BITS['A']+BITS['B'])
print('bits',len(allbits))
out={}
for b in allbits:
    pins=[]
    for a in chk[b]:
        p=polys[a]
        # look for pattern: coeff on (b,) is -c*HUGE ; monomial (b,F) coeff c
        m2=[(k,c) for k,c in p.items() if len(k)==2 and b in k]
        m1=[(k,c) for k,c in p.items() if k==(b,)]
        if len(m2)==1 and len(m1)==1 and len(p)<=3:
            k,c=m2[0]; F=k[0] if k[1]==b else k[1]
            if F==b: continue
            cc=m1[0][1]
            if cc % c: continue
            H=-cc//c
            pins.append((a,F,H,c))
    out[b]=pins
cnt=defaultdict(int)
for b,p in out.items(): cnt[len(p)]+=1
print('pins-per-bit hist',dict(cnt))
bad=[b for b,p in out.items() if len(p)!=2]
print('bits without exactly 2 pins:',bad[:20], len(bad))
json.dump({str(b):[(F,H,c) for a,F,H,c in v] for b,v in out.items()}, open(os.path.join(HERE,'leafpins.json'),'w'))
for b in allbits[:5]:
    print(b, out[b])

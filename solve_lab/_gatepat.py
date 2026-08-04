import heal_harness as H, re, json
from collections import Counter
pat=Counter()
ops=Counter()
for t,rhs,vids in H.gates:
    # classify
    r=rhs
    has_p = any(c.isdigit() for c in r)
    n_mul=r.count('*'); n_plus=r.count('+'); n_minus=r.count('-')
    key=(n_mul,n_plus,n_minus, len(vids))
    pat[key]+=1
# show top patterns
for k,c in pat.most_common(20):
    print(f"  (mul={k[0]},plus={k[1]},minus={k[2]},nvids={k[3]}): {c}")
# show a few example rhs per pattern
print("\nexamples:")
seen=set()
for t,rhs,vids in H.gates:
    n_mul=rhs.count('*'); n_plus=rhs.count('+'); n_minus=rhs.count('-')
    key=(n_mul,n_plus,n_minus,len(vids))
    if key not in seen and len(seen)<15:
        seen.add(key)
        print(f"  {key}: {rhs[:60]}")

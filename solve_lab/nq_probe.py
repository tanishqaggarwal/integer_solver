import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
# fill non-free from file too (gates get overwritten by forward)
for k,vv in d.items(): H.val[k]=vv
H.forward()
F=H.fails()
print(f"sy_regime11 fails ({len(F)}): {F}")
print(f"x_9062={H.val[9062]}, x_4287={H.val[4287]}, x_15298={H.val[15298]}")
# For each failing eq, show its residual and granularity (is residual a multiple of p?)
ns={'v':H.val,'__builtins__':{}}
for i in F:
    r=eval(H.eqcode[i],ns)
    gp = "p-MULTIPLE" if r%p==0 else f"sub-p (r%p has {len(str(r%p))} digits)"
    print(f"  eq{i}: residual {gp}, |r|~{len(str(abs(r)))} digits")

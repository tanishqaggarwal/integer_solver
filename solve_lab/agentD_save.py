#!/usr/bin/env python3
"""Save a specified config's full solution to JSON for independent checker verification."""
import sys, json
import agentD_harness as H
C1,C2=H.C1,H.C2
mode=sys.argv[1] if len(sys.argv)>1 else '10'
out=sys.argv[2] if len(sys.argv)>2 else 'agentD_test.json'
CONST={30213:C2, 22162:C1, 24468:C1, 18956:C2}
a7=24601; a34=2081
if mode=='11': ov={a7:1,a34:1,**CONST}
elif mode=='10': ov={a7:1,**CONST}
elif mode=='01': ov={a34:1,**CONST}
elif mode=='00': ov={**CONST}
elif mode=='00noC': ov={}
else: ov={a7:1,**CONST}
r=H.run_config(ov, want_val=True)
print("sat",r['satisfied'],"nfail",r['nfail'],"F",r['F'])
val=r['val']
json.dump({f"x_{i}":val[i] for i in range(H.NVARS)}, open(out,'w'))
print("saved",out)

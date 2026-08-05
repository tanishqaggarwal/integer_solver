import json
import heal_harness as H
p=H.p
# boolean/bit variables: the 256 pin selectors + bool_vars.json
pinrec=json.load(open('pinrec.json'))
bits=set(sel for i,sel,tgt,const,coef,handle in pinrec)
try:
    bv=json.load(open('atoms/bool_vars.json'))
    if isinstance(bv,dict): bits|=set(int(k) if not str(k).startswith('x') else int(k[2:]) for k in bv)
    elif isinstance(bv,list): bits|=set(int(x[2:]) if isinstance(x,str) else int(x) for x in bv)
except Exception as e: print("bool_vars:",e)
print(f"total bit/boolean vars: {len(bits)}")
# classify multiply gates
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gates.append((dd['t'],dd['rhs'],tuple(dd['vids'])))
import re
VAR=re.compile(r'x_(\d+)')
mult_gates=[]  # gates with a product
for t,rhs,vids in gates:
    # detect product: rhs contains '*' between two x_ vars
    if '*' in rhs:
        vs=[int(m) for m in VAR.findall(rhs)]
        mult_gates.append((t,rhs,vs))
print(f"total gates: {len(gates)}, multiply gates: {len(mult_gates)}")
# value*value: a product x_a*x_b where NEITHER a nor b is a bit AND neither is a constant
vv=0; bv_cnt=0; other=0
vv_gates=[]
for t,rhs,vs in mult_gates:
    # find the two factors of the product (crude: vars in rhs)
    # a gate rhs like "x_a * x_b" or "coef * (x_a * x_b)" etc.
    prodvars=[v for v in vs]
    nonbit=[v for v in prodvars if v not in bits]
    if len([v for v in prodvars if v not in bits])>=2:
        # could be value*value; check it's a genuine 2-var product
        vv+=1; vv_gates.append((t,rhs))
    elif any(v in bits for v in prodvars):
        bv_cnt+=1
    else: other+=1
print(f"value*value gates (genuine nonlinearity): {vv}")
print(f"bit*value gates (linear given bits): {bv_cnt}")
print("sample value*value gates:")
for t,rhs in vv_gates[:15]: print(f"  x_{t} = {rhs[:60]}")

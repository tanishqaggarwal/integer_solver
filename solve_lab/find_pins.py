import json
p=2**256-2**32-977
with open('atoms/poly_atoms.jsonl') as f:
    lines=f.readlines()
# pin pattern: sel*(x_i - BIGCONST) - k*handle
# = sel*x_i - sel*BIGCONST*... actually expanded: terms with a big constant coefficient
pins=[]
for i,l in enumerate(lines):
    a=json.loads(l)
    poly=a['poly']
    # look for a term that is [[selvar,var], -BIGCONST] i.e. product of sel*var with huge coeff? 
    # Actually repr form: sel * (x_i - CONST) - k*handle. Detect big const in repr.
    r=a['repr']
    # find constants in poly with |c| > p/4 (message-sized)
    bigc=[t for t in poly if not t[0] and abs(t[1])>p//4]  # constant terms
    bigcoef=[t for t in poly if t[0] and abs(t[1])>p//4]   # var terms with huge coeff
    if bigcoef and len(poly)<=5:
        pins.append((i,a))
print(f"candidate pin/message atoms (huge var-coeff, small): {len(pins)}")
for i,a in pins[:40]:
    print(f"atom{i} (n_eq={a['n_eq']}): {a['repr'][:130]}")

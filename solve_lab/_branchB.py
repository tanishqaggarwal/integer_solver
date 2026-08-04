import heal_harness as H, json, re
p=H.p
VAR=re.compile(r'x_(\d+)')
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f: gates.append(json.loads(line))
byout={g['t']:g for g in gates}
def defn(v):
    g=byout.get(v); return g['rhs'] if g else None
for v in [2099,19964,9062,4287,2081,21279,9118,8731]:
    r=defn(v)
    print(f"x_{v} <= {r[:150] if r else '(FREE INPUT)'}")
# trace x_2099 and x_19964 back a few levels
def expand(v,depth,maxd=4):
    ind='  '*depth; r=defn(v)
    if r is None: 
        tag=''
        print(ind+f"x_{v} [FREE]"); return
    if depth>=maxd: print(ind+f"x_{v} <= {r[:60]}..."); return
    print(ind+f"x_{v} <= {r[:90]}")
    for u in sorted(set(int(m) for m in VAR.findall(r))):
        expand(u,depth+1,maxd)
print("\n=== x_2099 tree ===")
expand(2099,0,maxd=4)
print("\n=== x_19964 tree ===")
expand(19964,0,maxd=4)

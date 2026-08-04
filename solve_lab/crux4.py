import pickle, json
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
cl=pickle.load(open(SCR+'/classify22.pkl','rb')); nonlin=set(cl['nonlin']); affine=set(cl['affine'])
cd=pickle.load(open(SCR+'/cert_nobool.pkl','rb')); cert=cd['cert']
reprs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        reprs[i]=json.loads(line).get('repr','')
print("Full cert (53 atoms), the 4 NONLINEAR ones:")
for a,mv in cert:
    if a in nonlin:
        print(f"  atom {a}: {reprs[a][:130]}")
print("\nAll cert atoms (id): ", sorted(a for a,_ in cert))
print("\nAffine cert atoms sample:")
for a,mv in cert[:20]:
    if a in affine:
        print(f"  atom {a}: {reprs[a][:90]}")

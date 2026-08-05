import pickle, json
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
cq=pickle.load(open(SCR+'/cert_quadrant.pkl','rb'))
reprs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f): reprs[i]=json.loads(line).get('repr','')
out={'quadrant':'x_2081=1,x_24601=1,MUX x_15298=1 (forward_construct overrides)',
     'certificate_atoms':[{'atom':a,'mult':str(mv),'repr':reprs[a]} for a,mv in sorted(cq['cert'])],
     'flip_candidates':[2081,24601],
     'gap_atoms':[20862,20864],'load_atoms':[18081,18084,29377,35321]}
json.dump(out, open('quadrant_infeasibility_certificate.json','w'), indent=1)
print("saved quadrant_infeasibility_certificate.json with", len(cq['cert']), "atoms")

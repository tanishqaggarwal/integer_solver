import heal_harness as H, sz_engine as E
p=H.p; RIP=set(E.RIP); CORE=set(E.CORE)
v013=H.loadd('best/new_instance_partial_39013.json')
v022=H.loadd('best_agentA_39022.json')
changed=[2498,2964,6083,11080,14623,14853,23238,24548,28246,31339,36462]  # core-frees (excl 4432,7068)
# endpoints: both at residues for x_7068,x_4432
def set_state(t):
    for v in H.freeinp: H.val[v]=v013.get(v,0)
    for w in changed:
        a=v013.get(w,0); b=v022.get(w,0)
        H.val[w]=a+ (b-a)*t//100   # integer interpolation, t in 0..100
    H.forward()
    H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
print(" t%   total  core  ripple  other")
for t in [0,10,20,30,40,50,60,70,80,90,95,99,100]:
    set_state(t)
    F=set(H.fails())
    print(f" {t:3d}   {len(F):4d}  {len(F&CORE):4d}  {len(F&RIP):5d}  {len(F-CORE-RIP):5d}")
# endpoint check
set_state(100); print("\nt=100 (core-solved+residues) total",len(H.fails()),"== hybrid 16?",len(H.fails())==16)
set_state(0);   print("t=0   (39013)                 total",len(H.fails()),"== 20?",len(H.fails())==20)

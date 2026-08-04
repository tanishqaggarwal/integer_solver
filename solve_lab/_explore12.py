import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
def atomval(ai):
    a=atoms[ai]; s=0
    for vl,c in a['poly']:
        t=c
        for v in vl: t*=H.val[v]
        s+=t
    return s
d=H.loadd('best/new_instance_partial_39013.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d); H.forward()
x24908=H.val[24908]
d2=dict(d); d2[14853]=d[12186]; d2[16742]=x24908
setfree(d2); H.forward()
a25170=atomval(25170); a27902=atomval(27902)
a42851=atomval(42851); a43834=atomval(43834); a44270=atomval(44270)
print("atom25170 =",a25170%p)
print("atom27902 =",a27902%p)
print("--- verifier residuals after subtracting embedded handle multiple ---")
print("atom42851 - (-26)*atom27902 =", (a42851 - (-26)*a27902)%p, " (exact:", a42851-(-26)*a27902==0,")")
print("atom44270 - (-5)*atom27902  =", (a44270 - (-5)*a27902)%p,  " (exact:", a44270-(-5)*a27902==0,")")
print("atom43834 - (2)*atom25170   =", (a43834 - (2)*a25170)%p,   " (exact:", a43834-(2)*a25170==0,")")

import heal_harness as H, re, json, pickle
p=H.p
# Parse gate rhs into (op, operands) where operand = ('v',id) or ('c',const)
def parse_rhs(rhs):
    rhs=rhs.strip()
    for op,sym in [('mul','*'),('add','+'),('sub','-')]:
        # split on top-level sym (no parens in gate rhs)
        if sym in rhs:
            parts=rhs.split(sym)
            if len(parts)==2:
                def opnd(s):
                    s=s.strip()
                    m=re.fullmatch(r'x_(\d+)',s)
                    if m: return ('v',int(m.group(1)))
                    return ('c',int(s))
                return (op,opnd(parts[0]),opnd(parts[1]))
    m=re.fullmatch(r'x_(\d+)',rhs)
    if m: return ('id',('v',int(m.group(1))))
    return ('id',('c',int(rhs)))
gate_parsed={}
for t,rhs,vids in H.gates:
    gate_parsed[t]=parse_rhs(rhs)
# reverse-mode AD: return grad over all vars for output var 'out'
def grad(out, val):
    a=[0]*H.NVARS
    a[out]=1
    for t in reversed(H.order):
        at=a[t]
        if at==0: continue
        pr=gate_parsed[t]
        op=pr[0]
        if op=='id':
            k=pr[1]
            if k[0]=='v': a[k[1]]=(a[k[1]]+at)%p
        elif op=='add':
            for k in (pr[1],pr[2]):
                if k[0]=='v': a[k[1]]=(a[k[1]]+at)%p
        elif op=='sub':
            k1,k2=pr[1],pr[2]
            if k1[0]=='v': a[k1[1]]=(a[k1[1]]+at)%p
            if k2[0]=='v': a[k2[1]]=(a[k2[1]]-at)%p
        elif op=='mul':
            k1,k2=pr[1],pr[2]
            v1 = val[k1[1]] if k1[0]=='v' else k1[1]
            v2 = val[k2[1]] if k2[0]=='v' else k2[1]
            if k1[0]=='v': a[k1[1]]=(a[k1[1]]+at*v2)%p
            if k2[0]=='v': a[k2[1]]=(a[k2[1]]+at*v1)%p
    return a
d0=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d0.get(v,0)
H.forward()
val=list(H.val)
gS=grad(35389,val); gT=grad(6671,val)
frees=sorted(H.freeinp)
nzS=[(f,gS[f]) for f in frees if gS[f]%p!=0]
nzT=[(f,gT[f]) for f in frees if gT[f]%p!=0]
print(f"free inputs with dS/df != 0: {len(nzS)}")
print(f"free inputs with dT/df != 0: {len(nzT)}")
print("frees moving S:", sorted(f for f,_ in nzS))
pickle.dump({'gate_parsed':gate_parsed,'gS':gS,'gT':gT,'val':val},
            open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/adcache.pkl','wb'))

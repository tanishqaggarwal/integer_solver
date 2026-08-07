"""CRT pairing + z-factor plumbing + beta->root permutation inversion."""
import pickle, math, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
import umodel as U
WZ=pickle.load(open('/home/user/integer_solver/solve_lab/agentU_work/w_z.pkl','rb'))
# z plumbing keyed by pin wire
ZINFO={}
for (sel,wire,z,m,fa,fb,ffa,ffb) in WZ['rows']:
    ZINFO[wire]=(z,m,fa,fb,ffa,ffb)

def crt(a1,m1,a2,m2):
    g=math.gcd(m1,m2)
    if (a2-a1)%g: return None
    l=m1//g*m2
    _,s,_=egcd(m1//g,m2//g)
    x=(a1 + (m1)*(((a2-a1)//g)*s % (m2//g)))%l
    return x
def egcd(a,b):
    if b==0: return (a,1,0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)

def pair_ok(sa,sb):
    wxa=[w for w in U.PIN[sa][0] if U.WCON[w][3]=='X'][0]
    wxb=[w for w in U.PIN[sb][0] if U.WCON[w][3]=='X'][0]
    ma,Ca=U.WCON[wxa][0],U.WCON[wxa][2]; mb,Cb=U.WCON[wxb][0],U.WCON[wxb][2]
    g=math.gcd(ma,mb)
    if (Ca-Cb)%g: return None
    W=crt(Ca,ma,Cb,mb)
    return W

def leaf_extras(sel, Wx, Wy):
    """z factor settings so the two pins hold over Z. Returns dict or None."""
    out={}
    for w in U.PIN[sel][0]:
        m,z,C,ax,e=U.WCON[w]
        val=Wx if ax=='X' else Wy
        d=val-C
        if d % m: return None
        q=d//m
        zi=ZINFO.get(w)
        if zi is None: return None
        zz,mm,fa,fb,ffa,ffb=zi
        if q==0:
            out[fa]=0; out[fb]=0
        elif ffa and ffb:
            out[fa]=q; out[fb]=1
        elif ffa:
            out[fa]=q; # fb pinned to p
            return None
        else:
            return None
    return out

# permutation from beta up to root: rootval[i] = betaval[chain(i)]
def chain_to_root(beta):
    idx=[0,1]
    n=beta
    while n in U.par:
        pnode=U.par[n]
        side='va' if U.tree[pnode][0]==n else 'vb'
        pm=U.perm[(pnode,side)]
        idx=[pm[i] for i in idx]   # parentval[i] = childval[pm[i]]
        n=pnode
    return idx     # rootval[i] = betaval[idx[i]]

def betaval_for(beta, tgt):
    idx=chain_to_root(beta)
    bv=[None,None]
    for i in (0,1): bv[idx[i]]=tgt[i]
    return tuple(bv)

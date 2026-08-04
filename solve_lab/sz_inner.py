"""Extract signed inner residual expressions by splitting top-level product factors."""
import re
VAR=re.compile(r'x_(\d+)')
def toplevel_factors(s):
    """split balanced expression s into top-level '*'-separated factors."""
    s=s.strip()
    facs=[]; depth=0; cur=''; i=0
    while i<len(s):
        c=s[i]
        if c=='(': depth+=1; cur+=c
        elif c==')': depth-=1; cur+=c
        elif c=='*' and depth==0:
            facs.append(cur); cur=''
        else: cur+=c
        i+=1
    facs.append(cur)
    return [f.strip() for f in facs]

def inner_expr(lhs):
    """Return (inner_python_expr, kind). kind: 'sq' if LHS=BASE*BASE, else 'lin' (product of const-ish * BASE)."""
    facs=toplevel_factors(lhs)
    # strip fully-paren wrapper factors that are constants (no x_)
    varfacs=[f for f in facs if VAR.search(f)]
    constfacs=[f for f in facs if not VAR.search(f)]
    if len(varfacs)==2 and varfacs[0]==varfacs[1]:
        base=varfacs[0]; kind='sq'
    elif len(varfacs)==1:
        base=varfacs[0]; kind='lin'
    else:
        # e.g. two DIFFERENT var factors (rare) -> product; treat as 'prod', inner = whole (fallback)
        base='('+')*('.join(varfacs)+')'; kind='prod'
    expr=VAR.sub(r'v[\1]',base)
    return expr,kind,len(constfacs)

if __name__=='__main__':
    import heal_harness as H, sz_engine as E
    from math import isqrt
    lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
    E.classify(); E.setup()
    ns={'v':H.val,'__builtins__':{}}
    ok=0; bad=[]
    for i in E.RIP:
        lhs=lines[i].rsplit('=',1)[0]
        expr,kind,ncf=inner_expr(lhs)
        code=compile(expr,'<i>','eval')
        Ev=eval(code,ns); LHSv=eval(H.eqcode[i],ns)
        # verify: sq -> LHSv == Ev^2 ; lin -> LHSv == const*Ev (const divides)
        if kind=='sq':
            good=(LHSv==Ev*Ev)
        else:
            good=(Ev!=0 and LHSv%Ev==0) or (Ev==0 and LHSv==0)
        ok+=good; 
        if not good: bad.append((i,kind))
        print(f"eq {i:6d} kind={kind} ncf={ncf} innerE={'0' if Ev==0 else str(Ev)[:20]+'..'} LHS==E^2/c?{good}")
    print(f"\nverified {ok}/{len(E.RIP)}  bad={bad}")

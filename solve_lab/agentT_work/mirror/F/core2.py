import sys,os
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from parse import parse_line,node_str,const_val
from parse2 import factors,flatten_sum_top
def core2(e):
    terms=[];flatten_sum_top(e,1,terms)
    tot=0; sset=set(); pw=set(); ref=None
    for sg,nd in terms:
        fs=[];factors(nd,fs); k=sg; nonc=[]
        for f in fs:
            cv=const_val(f)
            if cv is None: nonc.append(f)
            else: k*=cv
        ss=set(node_str(x) for x in nonc)
        if len(ss)!=1: return None,None,None
        sset|=ss; pw.add(len(nonc)); tot+=k; ref=nonc[0]
    if len(sset)!=1 or len(pw)!=1: return None,None,None
    return tot, list(pw)[0], ref

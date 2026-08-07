p='/home/user/integer_solver/solve_lab/agentO_work/simO.py'
s=open(p).read()
s=s.replace("    kw=dict(maxcore=maxcore,maxcorebits=maxcorebits)",
            "    if verbose: print('    system:',info,flush=True)\n    kw=dict(maxcore=maxcore,maxcorebits=maxcorebits)")
s=s.replace("        info['full_msg']=msg\n        sol,keep,blocked=maxsolvable",
            "        info['full_msg']=msg\n        if verbose: print('    full unsat:',msg,'-> greedy over',len(use),'rows',flush=True)\n        sol,keep,blocked=maxsolvable")
open(p,'w').write(s)
print('patched')

import json,os,numpy as np
d=json.load(open('xdata.json')); p=int(d['p'])
lad=[(int(a),int(b)) for a,b in d['ladder']]; T=(int(d['T'][0]),int(d['T'][1]))
M=(1<<64)-1
print('|S|=0 : T == O ?', False)
print('|S|=1 : T equals some 2^i*G ?', any(T==q for q in lad), '  (x-match:', any(T[0]==q[0] for q in lad),')')
k=np.memmap('tbl4s.bin',dtype=np.uint64,mode='r')
q=np.uint64(T[0]&M)
i=int(np.searchsorted(k,q))
print('b=0 case: TX low64 in table (|S|<=4 direct)?', i<len(k) and k[i]==q)

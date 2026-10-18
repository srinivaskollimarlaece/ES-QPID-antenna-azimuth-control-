import json, numpy as np
from tune_global import *
def r3(x):
    return float(f"{x:.5g}")
base=json.load(open('final_base.json')); eq=json.load(open('final_eq.json'))
allr={}; allr.update(base); allr.update(eq)
locked={}
for k,v in allr.items():
    tag,L=k.split('|'); L=float(L)
    vals={a:r3(b) for a,b in v['vals'].items()}
    os_=evalc(v['name'],L,vals,v['opts'])
    J=np.mean([o['sITAE']+W_OS*max(o['OS']-2,0)+W_U*o['su']+W_SS*max(o['ess']-0.02,0)+20*o['div'] for o in os_])
    locked[k]=dict(name=v['name'],vals=vals,opts=v['opts'],J=float(J),J_unrounded=v['J'])
    print(f"{tag:12s} L={L} J={v['J']:.4f} -> rounded {J:.4f}  {vals}")
json.dump(locked,open('locked.json','w'),indent=1)

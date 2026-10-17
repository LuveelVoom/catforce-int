import sys,re
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import decode
T={'I':lambda x,y:(x,y),'R90':lambda x,y:(-y,x),'R180':lambda x,y:(-x,-y),'R270':lambda x,y:(y,-x),
   'FX':lambda x,y:(x,-y),'FY':lambda x,y:(-x,y),'D1':lambda x,y:(y,x),'D2':lambda x,y:(-y,-x)}
def norm(c):
    mx=min(x for x,y in c);my=min(y for x,y in c);return frozenset((x-mx,y-my) for x,y in c)
def sym(c):
    c=norm(c);return frozenset(k for k,f in T.items() if norm({f(x,y) for x,y in c})==c)
def char(H):
    n=len(H)
    if n==8: return '.'
    if n==1: return '*'
    if n==2:
        if 'R180' in H: return 'x'
        return '+'              # one mirror (orthogonal or diagonal): 4 rotations
    if n==4:
        if 'FX' in H and 'FY' in H: return '/'   # needs 90-degree turn
        if 'D1' in H and 'D2' in H: return '|'
        if 'R90' in H: return '|'
    return '*'
out=[];rep=[]
ma=sys.argv[2]
for ln in open(sys.argv[1]):
    if ln.startswith('cat '):
        e=ln.split();c=decode(e[1]);H=sym(c)
        w=max(x for x,y in c)+1;h=max(y for x,y in c)+1
        e[2]=ma;e[3]=str(-(w//2));e[4]=str(-(h//2));e[5]=char(H)
        rep.append((e[1],len(H),e[5]))
        out.append(' '.join(e)+'\n')
    else: out.append(ln)
open(sys.argv[3],'w').writelines(out)
for r in rep: print(r)

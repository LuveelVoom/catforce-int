import sys,re,collections
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
R=tables("B3/S23-a5")
pat=decode("o$b2o$bo$o!"); P=pat|{(-x,y) for x,y in pat}
def ring(c): return {(x+dx,y+dy) for x,y in c for dx in (-1,0,1) for dy in (-1,0,1)}-set(c)
def run(cats,G=130):
    cur=set(P)
    for c in cats: cur|=c
    last=0; first=None; runmax=0; r=0
    for g in range(G):
        bad=any(not(c<=cur) or (ring(c)&cur) for c in cats)
        if bad:
            last=g; r+=1; runmax=max(runmax,r)
            if first is None: first=g
        else: r=0
        cur=step(cur,R)
    rest=cur-set().union(*cats) if cats else cur
    return first,last,runmax,frozenset(rest)
def per(rest,cats):
    # period of leftover (with catalysts present) up to 8
    cur=set(rest)|set().union(*cats); s0=frozenset(cur)
    for p in range(1,9):
        cur=step(cur,R)
        if frozenset(cur)==s0: return p
    return 0
alone=run([])[3]
print("alone:",len(alone))
sols=[]
for ln in open(sys.argv[1]):
    if ln.startswith('#'): continue
    cs=[]
    for m in re.finditer(r'cat line (\d+) \[[^\]]*\] as (\S+) top-left \((-?\d+),(-?\d+)\)',ln):
        c=decode(m.group(2)); x0=int(m.group(3)); y0=int(m.group(4))
        c={(x+x0,y+y0) for x,y in c}
        cs.append((int(m.group(1)),frozenset(c|{(-x,y) for x,y in c})))
    sols.append(cs)
C=collections.Counter(); ex={}
for i,cs in enumerate(sols):
    cats=[c for _,c in cs]
    f,l,rm,rest=run(cats)
    need=[run([c for j,c in enumerate(cats) if j!=k])[3]!=rest for k in range(len(cats))]
    # does a single catalyst alone already recover & give same result
    p=per(rest,cats)
    key=(len(rest),p,'both needed' if all(need) else 'one is filler' if any(need) else 'neither needed', 'short' if l-f<=4 else 'mid' if l-f<=25 else 'long')
    C[key]+=1; ex.setdefault(key,i+1)
print("leftover cells | period (0=not settled) | necessity | interaction span  -> count (example #)")
for k,v in sorted(C.items()): print(k,v,ex[k])

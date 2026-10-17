import sys; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
R=tables("B3/S23-a5")
pat=decode("o$b2o$bo$o!"); full=pat|{(-x,y) for x,y in pat}; pn=norm(full)
def ring(c): return {(x+dx,y+dy) for x,y in c for dx in (-1,0,1) for dy in (-1,0,1)}-set(c)
def analyse(fr):
    comps=comps8(fr)
    # the pattern: the component that equals the 8-cell pattern (may be merged with a touching catalyst)
    pc=[c for c in comps if norm(c)==pn]
    cats=[c for c in comps if norm(c)!=pn]
    info=[]
    cur=set(fr); first={}; missing={i:0 for i in range(len(cats))}; lastback={}
    for g in range(0,101):
        for i,c in enumerate(cats):
            ok = c<=cur and not (ring(c)&cur)
            if not ok:
                first.setdefault(i,g); missing[i]+=1
            else:
                if i in first and missing[i] and i not in lastback: pass
        cur=step(cur,R)
    final=cur
    intact=[c<=final and not (ring(c)&final) for c in cats]
    rest=final-set().union(*cats) if cats else final
    return dict(found_pattern=bool(pc), ncats=len(cats), first=[first.get(i) for i in range(len(cats))], missing=[missing[i] for i in range(len(cats))], intact=intact, ash=len(rest), pop=len(final), cats=cats)
if __name__=="__main__":
    fr=frames('in/results.rle')
    print("row frame | catalyst pieces | first touched (gen) | gens not intact (of 100) | intact at gen 100 | leftover cells besides catalysts")
    import collections
    summ=collections.Counter()
    for i,(r,f) in enumerate(fr):
        a=analyse(f)
        untouched=sum(1 for x in a['first'] if x is None)
        brief=sum(1 for m,x in zip(a['missing'],a['first']) if x is not None and m<=3)
        summ[(untouched,len(a['cats']))]+=1
        print(f"{r:3d} {i:5d} | {a['ncats']:2d} | {a['first']} | {a['missing']} | {['Y' if x else 'N' for x in a['intact']]} | {a['ash']}{'' if a['found_pattern'] else '  (pattern touching a catalyst at start)'}")
    print("frames by (catalyst pieces never touched in 100 gens, total pieces):", dict(summ))

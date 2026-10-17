import json, collections
import os; HENSEL = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hensel.json')))
def tables(r):
    birth=[0]*256; surv=[0]*256; mode=-1; i=0
    while i < len(r):
        ch=r[i]
        if ch=='/': i+=1; continue
        if ch in 'Bb': mode=0; i+=1; continue
        if ch in 'Ss': mode=1; i+=1; continue
        i+=1; neg=False
        if i<len(r) and r[i]=='-': neg=True; i+=1
        let=''
        while i<len(r) and r[i].isalpha() and r[i] not in 'BbSs': let+=r[i].lower(); i+=1
        dest = birth if mode==0 else surv
        for n,l,masks in HENSEL:
            if n!=ch: continue
            on = True if not let else ((l not in let) if neg else (l in let))
            if on:
                for m in masks: dest[m]=1
    return birth, surv
NB = {(-1,-1):0x80,(0,-1):0x40,(1,-1):0x20,(-1,0):0x10,(1,0):0x08,(-1,1):0x04,(0,1):0x02,(1,1):0x01}
def step(cells, R):
    birth, surv = R; m = collections.defaultdict(int)
    for (x,y) in cells:
        m[(x,y)] |= 0
        for (dx,dy),b in NB.items(): m[(x+dx,y+dy)] |= NB[(-dx,-dy)]
    return {c for c,mask in m.items() if (surv[mask] if c in cells else (mask and birth[mask]))}
def decode(t):
    cells=set(); x=y=0; n=0
    for ch in t:
        if ch.isdigit(): n=n*10+int(ch)
        elif ch in 'ob$':
            n=n or 1
            if ch=='o':
                for i in range(n): cells.add((x+i,y))
                x+=n
            elif ch=='b': x+=n
            else: y+=n; x=0
            n=0
        elif ch=='!': break
    return cells
def comps8(cells):
    seen=set(); res=[]
    for c in sorted(cells):
        if c in seen: continue
        st=[c]; comp=set()
        while st:
            p=st.pop()
            if p in seen: continue
            seen.add(p); comp.add(p)
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    q=(p[0]+dx,p[1]+dy)
                    if q in cells and q not in seen: st.append(q)
        res.append(comp)
    return res
def norm(c):
    mx=min(x for x,y in c); my=min(y for x,y in c); return frozenset((x-mx,y-my) for x,y in c)
def orients(c):
    out=set(); cur=set(c)
    for _ in range(4):
        cur={(-y,x) for x,y in cur}; out.add(norm(cur)); out.add(norm({(-x,y) for x,y in cur}))
    return out
def show(c):
    x0=min(x for x,y in c); y0=min(y for x,y in c); x1=max(x for x,y in c); y1=max(y for x,y in c)
    return "\n".join("".join("O" if (x,y) in c else "." for x in range(x0,x1+1)) for y in range(y0,y1+1))
def groups(vals, gap):
    vals=sorted(set(vals)); out=[]; cur=[vals[0]]
    for v in vals[1:]:
        if v-cur[-1]>gap: out.append((cur[0],cur[-1])); cur=[v]
        else: cur.append(v)
    out.append((cur[0],cur[-1])); return out
def frames(fn):
    body="".join(l for l in open(fn).read().split("\n") if not l.startswith("#") and not l.startswith("x ="))
    cells=decode(body); out=[]
    for ri,(r0,r1) in enumerate(groups([y for x,y in cells], 9)):
        rc={(x,y) for x,y in cells if r0<=y<=r1}
        for (c0,c1) in groups([x for x,y in rc], 9):
            out.append((ri,{(x-c0,y-r0) for x,y in rc if c0<=x<=c1}))
    return out

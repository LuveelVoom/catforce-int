import golly as g

def kth_bit(i, k):
    return (i >> k) & 1

s = ""

for i in range(512):
    for x in range(3):
        for y in range(3):
            g.setcell(x, y, kth_bit(i, 3 * x + y))
    g.run(1)
    s += "%i, " % g.getcell(1, 1)

g.setclipstr(s)

# Rebuilds the generator page: python3 generator/src/build.py (from any folder).
# Writes generator/catforce_in_generator.html and the copy in CatForce/.
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
old=open('old_generator_v8.html').read().split('\n')
core='\n'.join(old[356:357]+old[378:444]+old[667:854])
a='''    for (const m of members.filter(m => m.role === 'bait')) {'''
core=core.replace(a,'''    const alsoCells = [].concat(...members.filter(m => m.role === 'also').map(m => m.cells));
    if (alsoCells.length) { const ab = bbox(alsoCells); setup.also = {rle: encodeRle(alsoCells.map(c => [c[0] - ab.x0, c[1] - ab.y0])), x: ab.x0 + ox, y: ab.y0 + oy}; }
'''+a)
page=open('a_head.html').read()+'\n<script>\n/* ---------- rule tables, RLE and simulation helpers ---------- */\n'+core+'\n\n'+open('c_ui.js').read()+'</script>\n'
for out in ('../catforce_in_generator.html', '../../CatForce/catforce_in_generator.html'):
    open(out,'w').write(page)

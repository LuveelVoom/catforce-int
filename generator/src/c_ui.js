/* ================= page logic ================= */
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const esc = s => String(s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const isInt = v => /^-?\d+$/.test(String(v).trim());
const tr = v => String(v ?? '').trim();
function getPath(o, p) { return p.split('.').reduce((a, k) => a == null ? a : a[k], o); }
function setPath(o, p, v) { const ks = p.split('.'); const last = ks.pop(); ks.reduce((a, k) => a[k], o)[last] = v; }

/* ---------- symmetry of the whole search (matches LifeAPI.h transforms, origin at the torus centre) ---------- */
const SYMS = [
  ['C1', 'None'],
  ['D2|', 'Mirror left–right, axis through x = 0'],
  ['D2|even', 'Mirror left–right, axis between x = −1 and 0'],
  ['D2-', 'Mirror top–bottom, axis through y = 0'],
  ['D2-even', 'Mirror top–bottom, axis between y = −1 and 0'],
  ['D2\\', 'Mirror on the diagonal y = x'],
  ['D2/', 'Mirror on the diagonal y = −x'],
  ['C2', 'Rotate 180° about cell (0,0)'],
  ['C2even', 'Rotate 180° about the corner at (−½,−½)'],
  ['C2|even', 'Rotate 180°, even width'],
  ['C2-even', 'Rotate 180°, even height'],
  ['C4', 'Rotate 90° about cell (0,0)'],
  ['C4even', 'Rotate 90° about the corner at (−½,−½)'],
  ['D4+', 'Both orthogonal mirrors, through cell (0,0)'],
  ['D4+even', 'Both orthogonal mirrors, even width and height'],
  ['D4+|even', 'Both orthogonal mirrors, even width'],
  ['D4+-even', 'Both orthogonal mirrors, even height'],
  ['D4x', 'Both diagonal mirrors'],
  ['D4xeven', 'Both diagonal mirrors, even'],
  ['D8', 'All 8 symmetries, about cell (0,0)'],
  ['D8even', 'All 8 symmetries, about the corner at (−½,−½)']
];
const TF = {
  RX: ([x, y]) => [x, -y], RXe: ([x, y]) => [x, -1 - y], RY: ([x, y]) => [-x, y], RYe: ([x, y]) => [-1 - x, y],
  R180: ([x, y]) => [-x, -y], R180e: ([x, y]) => [-1 - x, -1 - y], R180v: ([x, y]) => [-x, -1 - y], R180h: ([x, y]) => [-1 - x, -y],
  D: ([x, y]) => [y, x], AD: ([x, y]) => [-y, -x], ADe: ([x, y]) => [-1 - y, -1 - x],
  R90: ([x, y]) => [-y, x], R90e: ([x, y]) => [-1 - y, x]
};
const CHAIN = {
  'C1': [], 'D2|': ['RY'], 'D2|even': ['RYe'], 'D2-': ['RX'], 'D2-even': ['RXe'], 'D2\\': ['D'], 'D2/': ['AD'],
  'C2': ['R180'], 'C2even': ['R180e'], 'C2|even': ['R180h'], 'C2-even': ['R180v'],
  'C4': ['R90', 'R180'], 'C4even': ['R90e', 'R180e'],
  'D4+': ['RX', 'RY'], 'D4+even': ['RXe', 'RYe'], 'D4+|even': ['RYe', 'RX'], 'D4+-even': ['RXe', 'RY'],
  'D4x': ['D', 'AD'], 'D4xeven': ['D', 'ADe'], 'D8': ['R90', 'R180', 'D'], 'D8even': ['R90e', 'R180e', 'D']
};
const wrap64 = v => ((v + 32) % 64 + 64) % 64 - 32;
function symExpand(cells, sym) {
  let T = cells.map(c => [c[0], c[1]]);
  for (const h of CHAIN[sym] || []) T = T.concat(T.map(TF[h]));
  const seen = new Set(), out = [];
  for (const [x, y] of T) { const k = wrap64(x) + ',' + wrap64(y); if (!seen.has(k)) { seen.add(k); out.push([wrap64(x), wrap64(y)]); } }
  return out;
}

/* ---------- catalyst orientations ---------- */
// The 8 orientations: I, R90, R180, R270, flip y, flip x, transpose, anti-transpose
const OR = [p => [p[0], p[1]], p => [-p[1], p[0]], p => [-p[0], -p[1]], p => [p[1], -p[0]],
            p => [p[0], -p[1]], p => [-p[0], p[1]], p => [p[1], p[0]], p => [-p[1], -p[0]]];
// transforms CatForce applies for each catalyst symmetry character (CharToTransforms in CatForce.cpp)
const CHAR_TF = {'.': [0], '|': [0, 5], '-': [0, 4], '\\': [0, 6], '/': [0, 7], '+': [0, 1, 2, 3], '@': [0, 1, 2, 3], 'x': [0, 1, 4, 6], '*': [0, 1, 2, 3, 4, 5, 6, 7]};
const CHAR_NAMES = {'.': 'as drawn', '|': 'mirror left–right', '-': 'mirror top–bottom', '\\': 'diagonal mirror', '/': 'other diagonal mirror', '+': '4 rotations', 'x': 'rotation + mirrors', '*': 'all 8'};
function normKey(cells) {
  let mx = Infinity, my = Infinity;
  for (const [x, y] of cells) { if (x < mx) mx = x; if (y < my) my = y; }
  return cells.map(([x, y]) => (x - mx) + ',' + (y - my)).sort().join(';');
}
function symGroup(cells) { const k = normKey(cells); return OR.map((f, i) => normKey(cells.map(f)) === k ? i : -1).filter(i => i >= 0); }
function autoChar(H) {
  const n = H.length;
  if (n === 8) return '.';
  if (n === 1) return '*';
  if (n === 2) return H.includes(2) ? 'x' : '+';
  if (n === 4) return H.includes(4) && H.includes(5) ? '/' : '|';
  return '*';
}
function orientCount(cells, ch) { return new Set((CHAR_TF[ch] || [0]).map(i => normKey(cells.map(OR[i])))).size; }

const LIB = [
  ['Block', '2o$2o!'], ['Beehive', 'b2o$o2bo$b2o!'], ['Loaf', 'b2o$o2bo$bobo$2bo!'], ['Boat', '2o$obo$bo!'],
  ['Ship', '2o$obo$b2o!'], ['Tub', 'bo$obo$bo!'], ['Pond', 'b2o$o2bo$o2bo$b2o!'], ['Long boat', '2o$obo$bobo$2bo!'],
  ['Eater 1', '2o$obo$2bo$2b2o!'], ['Snake', '2obo$ob2o!'], ['Barge', 'bo$obo$bobo$2bo!'],
  ['Aircraft carrier', '2o$o2bo$2b2o!'], ['Long ship', '2o$obo$bobo$2b2o!'], ['Mango', 'b2o$o2bo$bo2bo$2b2o!'],
  ['Blinker', '3o!'], ['Toad', 'b3o$3o!'], ['Beacon', '2o$o$3bo$2b2o!']
];
const GOAL_CHIPS = [['Glider', 'bo$2bo$3o!'], ['Block', '2o$2o!'], ['Herschel', 'o$3o$obo$2bo!']];

/* ---------- state ---------- */
function defaultState() {
  return {
    rule: 'B3/S23', paste: '', sym: 'C1', patRaw: '', patDx: '0', patDy: '0', baits: [],
    cats: [], goals: [], maxJunk: '',
    maxGen: '100', numCat: '2', stable: '15', areas: [],
    startGen: '', lastGen: '', numTransparent: '', stopAfter: '', maxCatSize: '', fitW: '', fitH: '', alsoRle: '', alsoX: '0', alsoY: '0',
    output: 'results.rle', outputList: '', fullReport: '', threads: '', maxTime: '', maxResults: '', part: '', outputGap: '',
    includes: [], keepDup: false, keepInteract: false, verbose: false, splitGap: '6', splitMerge: false
  };
}
function newCat() {
  return {name: '', rle: '', maxActive: '20', sym: 'auto', autoOff: true, dx: '0', dy: '0',
    forbidden: [], req: {rle: '', x: '0', y: '0'}, anti: {rle: '', x: '0', y: '0'}, locus: {rle: '', x: '0', y: '0'},
    transparent: false, mustinclude: false, checkRecovery: false, sacrificial: false, isBlinker: false, period: ''};
}
const newBait = () => ({rle: '', x: '0', y: '0', survive: '', absent: ''});
const newArea = mode => ({mode, x: '', y: '', w: '', h: ''});
const newGoal = type => ({type, gen: '', rle: '', dx: '0', dy: '0'});

const SIZE_FIELDS = [
  ['maxGen', 'Generations to simulate', 'max-gen', '250', 'How long each candidate is run.'],
  ['numCat', 'Catalysts per result', 'num-catalyst', '2', 'Exactly this many (1–5), plus their symmetric copies. Each extra catalyst multiplies the search time.'],
  ['stable', 'Recovery time', 'stable-interval', '15', 'Generations every catalyst must stay intact at the end to count as recovered.']
];
const LIMIT_FIELDS = [
  ['startGen', 'First contact no earlier than', 'start-gen', '1', 'Generation of the first catalyst contact.'],
  ['lastGen', 'First contact no later than', 'last-gen', 'max-gen − 1', ''],
  ['numTransparent', 'Transparent catalysts allowed', 'num-transparent', '100', 'Counts catalysts marked transparent.'],
  ['stopAfter', 'Give up after a catalyst breaks', 'stop-after-cats-destroyed', 'never', 'Generations to keep checking goals after a catalyst fails.'],
  ['maxCatSize', 'Results per output row', 'max-category-size', 'no limit', 'Further results with the same outcome are dropped.']
];
const RUN_FIELDS = [
  ['output', 'Results file', 'output', 'results.rle', 'Open it in Golly with the same rule.', false],
  ['outputList', 'Text list of results', 'output-list', 'none', 'Each result as catalyst names and positions.', false],
  ['fullReport', 'Unfiltered results file', 'full-report', 'none', 'Every result, ignoring the goals in step 4.', false],
  ['threads', 'CPU threads', 'threads', 'all cores', 'Same results, only faster.', true],
  ['maxTime', 'Time limit', 'max-time', 'none', 'Stops cleanly and saves: 90, 15m, 2h, 1d.', false],
  ['maxResults', 'Stop after this many results', 'max-results', 'no limit', '', true],
  ['part', 'Run one share of the search', 'part', 'whole search', 'K/N, e.g. 2/4 on the second of four computers.', false],
  ['outputGap', 'Gap between results', 'output-gap', '10', 'Empty cells between results in the RLE.', true]
];
function fieldHTML([k, label, kw, ph, help, numeric]) {
  return `<label class="fld"><span class="lbl">${esc(label)} <code>${esc(kw)}</code></span>
    <input type="text" id="f-${k}" data-k="${k}" placeholder="${esc(ph)}" spellcheck="false"${numeric !== false ? ' inputmode="numeric"' : ''}>
    ${help ? `<small>${esc(help)}</small>` : ''}</label>`;
}

let S = defaultState();
let split = {objs: [], setups: [], pick: 0};
let catCache = new Map();

function ruleOk() { return ruleTables(S.rule); }
function ruleIsotropic() { const rc = checkRule(S.rule); return !(rc.ok && /not isotropic/.test(rc.note || '')); }

function catInfo(rle) {
  const key = S.rule + '|' + rle;
  if (catCache.has(key)) return catCache.get(key);
  const cl = cleanRle(rle), d = decodeRle(cl.rle);
  const info = {error: cl.error, cells: d.cells, w: d.w, h: d.h, n: d.cells.length};
  if (d.cells.length) {
    const H = symGroup(d.cells);
    info.full = 8 / H.length; info.auto = autoChar(H);
    const R = ruleOk();
    try { info.cls = R.ok && d.cells.length < 2000 ? classifyObject(new Set(d.cells.map(c => c.join(','))), R, 64) : null; } catch (e) { info.cls = null; }
  }
  catCache.set(key, info);
  return info;
}
function catResolved(c) {
  const info = catInfo(cleanRle(c.rle).rle);
  const sym = c.sym === 'auto' ? (ruleIsotropic() ? info.auto || '.' : '.') : c.sym;
  const dx = c.autoOff ? String(-(info.w >> 1)) : tr(c.dx), dy = c.autoOff ? String(-(info.h >> 1)) : tr(c.dy);
  return {info, sym, dx, dy, tries: info.n ? orientCount(info.cells, sym) : 0};
}
function catLabel(c, i) { return c.name || 'Catalyst ' + (i + 1); }

/* ---------- building the file ---------- */
function build(S, forSetup) {
  const L = [], issues = [];
  const err = m => issues.push({t: 'err', m}), warn = m => issues.push({t: 'warn', m}), tip = m => issues.push({t: 'tip', m});
  const num = (kw, v, min = 0) => {
    v = tr(v); if (!v) return;
    if (!isInt(v) || +v < min) { err(kw + ': "' + v + '" must be a whole number' + (min > 0 ? ' of at least ' + min : '')); return; }
    L.push(kw + ' ' + v);
  };
  const rleOf = (raw, what) => { const c = cleanRle(raw); if (c.error) err(what + ': ' + c.error); return c.rle; };
  const fname = v => tr(v).replace(/\s+/g, '_');

  const rule = String(S.rule || '').replace(/\s+/g, '');
  const rc = checkRule(rule);
  if (!rc.ok) err('Rule: ' + rc.error);
  if (rule) L.push('rule ' + rule);

  num('max-gen', S.maxGen, 1); num('num-catalyst', S.numCat, 1); num('stable-interval', S.stable);
  num('start-gen', S.startGen); num('last-gen', S.lastGen); num('num-transparent', S.numTransparent);
  num('stop-after-cats-destroyed', S.stopAfter); num('max-category-size', S.maxCatSize);
  const g = (v, d) => tr(v) === '' || !isInt(v) ? d : +v;
  const maxGen = g(S.maxGen, 250), startGen = g(S.startGen, 1), lastGen = g(S.lastGen, maxGen - 1), nc = g(S.numCat, 2), stable = g(S.stable, 15);
  if (startGen > lastGen) warn('The first contact window is empty (start-gen ' + startGen + ' is after last-gen ' + lastGen + ').');
  if (nc > 5) err('num-catalyst: CatForce places at most 5 catalysts.');
  else if (nc >= 4) tip(nc + ' catalysts per result can take days. Try a smaller number first.');
  if (stable >= maxGen) warn('Recovery time (' + stable + ') must be shorter than the generations simulated (' + maxGen + ').');

  if (S.sym && S.sym !== 'C1') {
    if (!CHAIN[S.sym]) err('Unknown symmetry "' + S.sym + '".');
    L.push('symmetry ' + S.sym);
  }
  S.areas.forEach((r, i) => {
    const v = [r.x, r.y, r.w, r.h].map(tr);
    if (v.every(t => t === '')) return;
    const n = (r.mode === 'exclude' ? 'Excluded rectangle ' : 'Search rectangle ') + (i + 1);
    if (v.some(t => !isInt(t))) { err(n + ' needs four whole numbers: x, y, width, height.'); return; }
    if (+v[2] <= 0 || +v[3] <= 0) { err(n + ': width and height must be positive.'); return; }
    if (r.mode !== 'exclude' && (+v[0] < -32 || +v[1] < -32 || +v[0] + +v[2] > 32 || +v[1] + +v[3] > 32)) warn(n + ' reaches past the edge of the torus (−32 to 31).');
    L.push((r.mode === 'exclude' ? 'exclude-area ' : 'search-area ') + v.join(' '));
  });

  const pat = rleOf(S.patRaw, 'Active pattern');
  if (!pat || pat === '!') err('There is no active pattern yet (step 2).');
  else if (!isInt(S.patDx) || !isInt(S.patDy)) err('The active pattern position must be whole numbers.');
  else L.push('pat ' + pat + ' ' + tr(S.patDx) + ' ' + tr(S.patDy));

  if (S.baits.length > 8) err('CatForce accepts at most 8 baits.');
  S.baits.forEach((b, i) => {
    const n = 'Bait ' + (i + 1), rle = rleOf(b.rle, n);
    if (!rle) { err(n + ' has no pattern.'); return; }
    if (!isInt(b.x) || !isInt(b.y)) { err(n + ': x and y must be whole numbers.'); return; }
    let line = 'bait ' + rle + ' ' + tr(b.x) + ' ' + tr(b.y);
    for (const k of ['survive', 'absent']) {
      const t = tr(b[k]); if (!t) continue;
      if (!/^\d+$/.test(t)) { err(n + ': ' + k + ' must be a whole number.'); return; }
      line += ' ' + k + ' ' + t;
    }
    L.push(line);
  });

  const also = rleOf(S.alsoRle, 'Always present');
  if (also) {
    if (!isInt(S.alsoX) || !isInt(S.alsoY)) err('Always present: x and y must be whole numbers.');
    else L.push('also-required ' + also + ' ' + tr(S.alsoX) + ' ' + tr(S.alsoY));
  }
  if (tr(S.fitW) || tr(S.fitH)) {
    if (!isInt(S.fitW) || !isInt(S.fitH)) err('Fit in width and height needs both numbers.');
    else L.push('fit-in-width-height ' + tr(S.fitW) + ' ' + tr(S.fitH));
  }

  // goals
  if (S.goals.length) L.push('', '# goal');
  S.goals.forEach((f, i) => {
    const n = 'Goal ' + (i + 1), rle = rleOf(f.rle, n), gen = tr(f.gen) || defaultGoalGen();
    if (!rle) { err(n + ' has no pattern.'); return; }
    if (!/^\d+(-\d+)?$/.test(gen)) { err(n + ': generation must be a number such as 60 or a range such as 20-150.'); return; }
    if (f.type === 'match') { L.push('match ' + gen + ' ' + rle); return; }
    if (!isInt(f.dx) || !isInt(f.dy)) { err(n + ': x and y must be whole numbers.'); return; }
    L.push(f.type + ' ' + gen + ' ' + rle + ' ' + tr(f.dx) + ' ' + tr(f.dy));
  });
  num('max-junk', S.maxJunk);
  if (!S.goals.length && !S.baits.length) tip('No goal (step 4): every reaction where the catalysts recover is reported, including ones where the pattern just dies.');

  // output and running
  if (tr(S.output)) L.push('', 'output ' + fname(S.output));
  if (tr(S.outputList)) L.push('output-list ' + fname(S.outputList));
  if (tr(S.fullReport)) L.push('full-report ' + fname(S.fullReport));
  num('threads', S.threads);
  if (tr(S.maxTime)) { if (!/^\d+[smhd]?$/i.test(tr(S.maxTime))) err('Time limit must look like 90, 15m, 2h or 1d.'); else L.push('max-time ' + tr(S.maxTime)); }
  num('max-results', S.maxResults, 1);
  if (tr(S.part)) {
    const m = tr(S.part).match(/^(\d+)\/(\d+)$/);
    if (!m || +m[1] < 1 || +m[1] > +m[2]) err('Share of the search must look like 2/4 (K from 1 to N).');
    else L.push('part ' + tr(S.part));
  }
  num('output-gap', S.outputGap);
  if (S.keepDup) L.push('keep-duplicates yes');
  if (S.keepInteract) L.push('keep-interacting yes');
  if (S.verbose) L.push('verbose yes');
  const includes = S.includes.map(tr).filter(Boolean);
  includes.forEach(x => L.push('include ' + x));

  // catalysts
  if (!S.cats.length && !includes.length) err('Add at least one catalyst (step 3).');
  if (S.cats.length) L.push('', '# catalysts');
  let fewOrient = 0, longMissing = 0;
  S.cats.forEach((c, i) => {
    const n = catLabel(c, i), rle = rleOf(c.rle, n);
    if (!rle || rle === '!') { err(n + ' has no pattern.'); return; }
    const rs = catResolved(c), info = rs.info;
    if (!isInt(c.maxActive) || +c.maxActive < 0) { err(n + ': "may be missing for" must be a whole number.'); return; }
    if (!isInt(rs.dx) || !isInt(rs.dy)) { err(n + ': offset must be whole numbers.'); return; }
    if (!CHAR_TF[rs.sym]) { err(n + ': unknown orientation code "' + rs.sym + '".'); return; }
    let line = 'cat ' + rle + ' ' + tr(c.maxActive) + ' ' + rs.dx + ' ' + rs.dy + ' ' + rs.sym;
    const xy = (kw, o) => {
      const r2 = rleOf(o.rle, n + ' ' + kw); if (!r2) return '';
      if (!isInt(o.x) || !isInt(o.y)) { err(n + ' ' + kw + ': x and y must be whole numbers.'); return ''; }
      return ' ' + kw + ' ' + r2 + ' ' + tr(o.x) + ' ' + tr(o.y);
    };
    for (const f of c.forbidden) line += xy('forbidden', f);
    line += xy('required', c.req) + xy('antirequired', c.anti) + xy('locus', c.locus);
    if (c.transparent) line += ' transparent';
    if (c.mustinclude) line += ' mustinclude';
    if (c.checkRecovery) line += ' check-recovery';
    if (c.sacrificial) line += ' sacrificial';
    if (c.isBlinker) line += ' isblinker';
    if (tr(c.period)) { if (!isInt(c.period) || +c.period < 1) err(n + ': period must be a positive whole number.'); else line += ' period ' + tr(c.period); }
    L.push(line);
    const cls = info.cls;
    if (cls && !c.sacrificial) {
      if (cls.type === 'osc' && !tr(c.period)) warn(n + ' oscillates with period ' + cls.period + '. Set its period (step 3) or it can never count as recovered.');
      else if (cls.type !== 'still' && cls.type !== 'osc') warn(n + ' is not stable in this rule (' + cls.label + '), so it can never recover.');
    }
    if (rs.tries < info.full) fewOrient++;
    if (+c.maxActive >= 50) longMissing++;
  });
  if (fewOrient) warn(fewOrient + ' catalyst' + (fewOrient > 1 ? 's are' : ' is') + ' only tried in some orientations. Use "Try every orientation" in step 3 unless that is deliberate.');
  if (longMissing) tip(longMissing + ' catalyst' + (longMissing > 1 ? 's' : '') + ' may be missing for 50 or more generations, so a destroyed catalyst that happens to reform still counts. 10 to 30 is typical.');
  return {text: L.join('\n').replace(/^\n+/, '') + '\n', issues};
}
function defaultGoalGen() { const mg = isInt(S.maxGen) ? +S.maxGen : 250; return '1-' + mg; }

/* ---------- reading an existing file ---------- */
function parseIn(text) {
  const st = defaultState(), ignored = [];
  Object.assign(st, {maxGen: '', numCat: '', stable: '', output: '', patDx: '0', patDy: '0', rule: ''});
  const simple = {'max-gen': 'maxGen', 'start-gen': 'startGen', 'last-gen': 'lastGen', 'num-catalyst': 'numCat', 'num-transparent': 'numTransparent',
    'stable-interval': 'stable', 'stop-after-cats-destroyed': 'stopAfter', 'max-category-size': 'maxCatSize', 'output-gap': 'outputGap',
    'threads': 'threads', 'max-results': 'maxResults', 'max-time': 'maxTime', 'max-junk': 'maxJunk', 'part': 'part', 'symmetry': 'sym'};
  const yes = v => !/^(no|false|0|off)$/i.test(v);
  for (const line of text.split(/\r?\n/)) {
    const t = line.trim().split(/\s+/), kw = t[0];
    if (!kw || !/^[a-z]/i.test(kw)) continue;
    if (t.length < 2) { ignored.push(line.trim()); continue; }
    if (simple[kw]) { st[simple[kw]] = t[1]; continue; }
    switch (kw) {
      case 'rule': st.rule = t.slice(1).join(''); break;
      case 'output': st.output = t.slice(1).join(' '); break;
      case 'output-list': st.outputList = t.slice(1).join(' '); break;
      case 'full-report': st.fullReport = t.slice(1).join(' '); break;
      case 'include': st.includes.push(t.slice(1).join(' ')); break;
      case 'keep-duplicates': st.keepDup = yes(t[1]); break;
      case 'keep-interacting': st.keepInteract = yes(t[1]); break;
      case 'verbose': st.verbose = yes(t[1]); break;
      case 'search-area': case 'exclude-area':
        if (t.length >= 5) st.areas.push({mode: kw === 'exclude-area' ? 'exclude' : 'include', x: t[1], y: t[2], w: t[3], h: t[4]}); break;
      case 'pat': st.patRaw = t[1]; if (t.length > 3) { st.patDx = t[2]; st.patDy = t[3]; } break;
      case 'fit-in-width-height': st.fitW = t[1]; st.fitH = t[2] || ''; break;
      case 'also-required': st.alsoRle = t[1]; st.alsoX = t[2] || '0'; st.alsoY = t[3] || '0'; break;
      case 'filter': case 'andfilter': case 'orfilter':
        st.goals.push({type: kw === 'filter' ? 'andfilter' : kw, gen: t[1], rle: t[2] || '', dx: t[3] || '0', dy: t[4] || '0'}); break;
      case 'match': st.goals.push({type: 'match', gen: t[1], rle: t[2] || '', dx: '0', dy: '0'}); break;
      case 'bait': {
        if (t.length < 4) { ignored.push(line.trim()); break; }
        const b = newBait(); b.rle = t[1]; b.x = t[2]; b.y = t[3];
        for (let k = 4; k < t.length; k++) { if (t[k] === 'survive') b.survive = t[++k] || ''; else if (t[k] === 'absent') b.absent = t[++k] || ''; }
        st.baits.push(b); break;
      }
      case 'cat': {
        if (t.length < 6) { ignored.push(line.trim()); break; }
        const c = newCat();
        c.rle = t[1]; c.maxActive = t[2]; c.dx = t[3]; c.dy = t[4]; c.sym = t[5][0] === '@' ? '+' : t[5][0]; c.autoOff = false;
        let i = 6;
        while (i < t.length) {
          const w = t[i];
          if (w === 'forbidden') { c.forbidden.push({rle: t[i + 1] || '', x: t[i + 2] || '0', y: t[i + 3] || '0'}); i += 4; }
          else if (w === 'required' || w === 'antirequired' || w === 'locus') {
            const o = w === 'required' ? c.req : w === 'antirequired' ? c.anti : c.locus;
            o.rle = t[i + 1] || ''; o.x = t[i + 2] || '0'; o.y = t[i + 3] || '0'; i += 4;
          }
          else if (w === 'period') { c.period = t[i + 1] || ''; i += 2; }
          else if (w === 'transparent') { c.transparent = true; i++; }
          else if (w === 'mustinclude') { c.mustinclude = true; i++; }
          else if (w === 'check-recovery') { c.checkRecovery = true; i++; }
          else if (w === 'sacrificial') { c.sacrificial = true; i++; }
          else if (w === 'isblinker') { c.isBlinker = true; i++; }
          else { ignored.push('catalyst option "' + w + '"'); i++; }
        }
        st.cats.push(c); break;
      }
      default: ignored.push(line.trim().slice(0, 60));
    }
  }
  return {st, ignored};
}

/* ---------- rendering ---------- */
function cssVar(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
function drawCells(cv, cells, color) {
  const ctx = cv.getContext('2d');
  ctx.clearRect(0, 0, cv.width, cv.height);
  if (!cells.length) return;
  const b = bbox(cells), sc = Math.max(2, Math.min(18, Math.floor((cv.width - 12) / Math.max(b.w, b.h))));
  const ox = Math.floor((cv.width - b.w * sc) / 2), oy = Math.floor((cv.height - b.h * sc) / 2);
  ctx.fillStyle = color || cssVar('--cell');
  for (const [x, y] of cells) ctx.fillRect(ox + (x - b.x0) * sc, oy + (y - b.y0) * sc, sc - (sc > 3 ? 1 : 0), sc - (sc > 3 ? 1 : 0));
}

function initStatic() {
  $('#f-sym').innerHTML = SYMS.map(([v, l]) => `<option value="${esc(v)}">${esc(l)}${v === 'C1' ? '' : '  (' + esc(v) + ')'}</option>`).join('');
  $('#sizeFields').innerHTML = SIZE_FIELDS.map(fieldHTML).join('');
  $('#limitFields').innerHTML = LIMIT_FIELDS.map(fieldHTML).join('');
  $('#runFields').innerHTML = RUN_FIELDS.map(fieldHTML).join('');
}
function syncStatic() {
  $$('[data-k]').forEach(el => {
    if (el.closest('.list')) return;
    const v = getPath(S, el.dataset.k);
    if (el.type === 'checkbox') el.checked = !!v;
    else if (el.tagName === 'SELECT') {
      if (v && ![...el.options].some(o => o.value === v)) el.add(new Option(v + ' (not recognised)', v));
      el.value = v || 'C1';
    } else el.value = v ?? '';
  });
  $('#f-includes').value = S.includes.join('\n');
}

function renderObjs() {
  const box = $('#objs');
  if (!split.objs.length) { box.innerHTML = ''; return; }
  box.innerHTML = split.objs.map((o, k) => `
    <div class="obj" data-k="${k}">
      <canvas width="128" height="128" aria-hidden="true"></canvas>
      <div><b>Object ${k + 1}</b> <span class="opt">${o.cells.length} cells, ${o.w}×${o.h}</span><br><span class="opt">${esc(o.info.label)}</span></div>
      <select class="role" aria-label="Role of object ${k + 1}">
        <option value="pattern"${o.role === 'pattern' ? ' selected' : ''}>Active pattern</option>
        <option value="bait"${o.role === 'bait' ? ' selected' : ''}>Bait</option>
        <option value="also"${o.role === 'also' ? ' selected' : ''}>Always present</option>
        <option value="skip"${o.role === 'skip' ? ' selected' : ''}>Leave out</option>
      </select>
    </div>`).join('');
  $$('.obj', box).forEach(row => {
    const o = split.objs[+row.dataset.k];
    drawCells($('canvas', row), o.cells, o.role === 'bait' ? cssVar('--bait') : o.role === 'also' ? cssVar('--warn') : o.role === 'skip' ? cssVar('--axis') : cssVar('--cell'));
  });
}
function renderSetups() {
  const box = $('#setups'), N = split.setups.length;
  if (!split.objs.length) { box.innerHTML = ''; return; }
  box.innerHTML = (N < 2 ? '' : `<p class="hint">These objects form ${N} separate groups, so they make ${N} independent searches. The form shows the one you pick; you can also export all of them.</p>` +
    split.setups.map((st, i) => `<label class="setup${i === split.pick ? ' on' : ''}"><input type="radio" name="setupPick" value="${i}"${i === split.pick ? ' checked' : ''}>
      <span><b>Search ${i + 1}</b> · ${st.n} object${st.n > 1 ? 's' : ''} · ${st.baits.length} bait${st.baits.length === 1 ? '' : 's'}${st.pat ? '' : ' · <span class="warn">no active pattern</span>'}</span></label>`).join('')) +
    `<div class="inline">
      <label class="fld"><span class="lbl">Same search if closer than</span><input type="text" id="f-splitGap" data-k="splitGap" inputmode="numeric" value="${esc(S.splitGap)}"><small>cells apart</small></label>
      <label class="check" style="align-self:center"><input type="checkbox" id="f-splitMerge" data-k="splitMerge"${S.splitMerge ? ' checked' : ''}><span>Keep everything as one search</span></label>
    </div>
    ` + (N < 2 ? '' : `<div class="btnrow"><button id="copyAll" class="small">Copy all ${N} files</button><button id="saveAll" class="small" ${dl ? '' : 'hidden'}>Save all ${N} files</button><span id="allMsg" class="msg"></span></div>`);
}
function renderPatSummary() {
  const box = $('#patSummary');
  const pc = cleanRle(S.patRaw), pd = decodeRle(pc.rle);
  const rows = [];
  rows.push(pd.cells.length
    ? `<div class="row"><b>Active</b><span>${pd.cells.length} cells, ${pd.w}×${pd.h}, top-left at (${esc(S.patDx || 0)}, ${esc(S.patDy || 0)})</span></div>`
    : `<div class="row"><b>Active</b><span class="err">none yet: paste a pattern above</span></div>`);
  S.baits.forEach((b, i) => {
    const d = decodeRle(cleanRle(b.rle).rle);
    rows.push(`<div class="row"><b>Bait ${i + 1}</b><span>${d.cells.length} cells at (${esc(b.x)}, ${esc(b.y)})${tr(b.survive) ? ', must survive ' + esc(b.survive) + ' gens' : ''}</span></div>`);
  });
  box.innerHTML = rows.join('');
}
function renderBaits() {
  $('#baits').innerHTML = S.baits.map((b, i) => `
    <div class="goal" data-i="${i}">
      <label class="fld grow"><span class="lbl">Bait ${i + 1} <code>bait</code></span>
        <textarea data-k="baits.${i}.rle" data-rle="1" rows="2" spellcheck="false">${esc(b.rle)}</textarea></label>
      <label class="fld"><span class="lbl">x</span><input type="text" data-k="baits.${i}.x" value="${esc(b.x)}" inputmode="numeric"></label>
      <label class="fld"><span class="lbl">y</span><input type="text" data-k="baits.${i}.y" value="${esc(b.y)}" inputmode="numeric"></label>
      <label class="fld"><span class="lbl">Survive <code>survive</code></span><input type="text" data-k="baits.${i}.survive" value="${esc(b.survive)}" inputmode="numeric" placeholder="recovery time"></label>
      <label class="fld"><span class="lbl">Gone at least <code>absent</code></span><input type="text" data-k="baits.${i}.absent" value="${esc(b.absent)}" inputmode="numeric" placeholder="1"></label>
      <button class="small danger" data-act="delBait" data-i="${i}">Remove</button>
    </div>`).join('') || '<div class="empty">No baits.</div>';
}

function renderLib() {
  const R = ruleOk();
  const have = new Set(S.cats.map(c => { const d = decodeRle(cleanRle(c.rle).rle); return d.cells.length ? canonKey(d.cells) : ''; }));
  $('#lib').innerHTML = LIB.map(([name, rle], i) => {
    const info = catInfo(rle), cls = info.cls;
    const ok = R.ok && cls && (cls.type === 'still' || cls.type === 'osc');
    const on = have.has(canonKey(info.cells));
    const sub = !ok ? 'unstable here' : cls.type === 'osc' ? 'p' + cls.period : '';
    return `<button class="chip${on ? ' on' : ''}" data-lib="${i}" ${ok && !on ? '' : 'disabled'} title="${esc(ok ? (on ? 'Already added' : 'Add ' + name) : name + ' is not a still life or oscillator in this rule')}">${esc(name)}${sub ? ' <span class="sub">' + esc(sub) + '</span>' : ''}</button>`;
  }).join('');
}
function canonKey(cells) { return OR.map(f => normKey(cells.map(f))).sort()[0]; }

function statusBadge(info) {
  const cls = info.cls;
  if (info.error) return `<span class="badge err">bad RLE</span>`;
  if (!info.n) return `<span class="badge warn">empty</span>`;
  if (!cls) return '';
  if (cls.type === 'still') return `<span class="badge ok">still life</span>`;
  if (cls.type === 'osc') return `<span class="badge ok">oscillator p${cls.period}</span>`;
  return `<span class="badge err">${esc(cls.label)}</span>`;
}
function subRow(path, o, label, i, j) {
  return `<div class="subrow"><textarea data-k="${path}.rle" data-rle="1" rows="1" spellcheck="false" aria-label="${esc(label)} pattern" placeholder="RLE">${esc(o.rle)}</textarea>
    <input type="text" data-k="${path}.x" value="${esc(o.x)}" aria-label="x" placeholder="x">
    <input type="text" data-k="${path}.y" value="${esc(o.y)}" aria-label="y" placeholder="y">
    ${j !== undefined ? `<button class="small danger" data-act="delForbid" data-i="${i}" data-j="${j}">Remove</button>` : '<span></span>'}</div>`;
}
function renderCats() {
  const openSet = new Set($$('#cats .cat details[open]').map(d => d.closest('.cat').dataset.i));
  let under = 0;
  $('#cats').innerHTML = S.cats.map((c, i) => {
    const rs = catResolved(c), info = rs.info, cls = info.cls;
    const unstable = info.n && cls && cls.type !== 'still' && cls.type !== 'osc' && !c.sacrificial;
    if (info.n && rs.tries < info.full) under++;
    const opts = [`<option value="auto"${c.sym === 'auto' ? ' selected' : ''}>Every orientation (${info.full || 1})</option>`]
      .concat(Object.keys(CHAR_NAMES).map(ch => {
        const k = info.n ? orientCount(info.cells, ch) : 1;
        return `<option value="${esc(ch)}"${c.sym === ch ? ' selected' : ''}>${esc(CHAR_NAMES[ch])}: ${k} of ${info.full || 1}  (${esc(ch)})</option>`;
      })).join('');
    let note = '';
    if (info.n && rs.tries < info.full) note += `<div class="cat-note warn">Only ${rs.tries} of its ${info.full} orientations are tried. <button class="small" data-act="autoSym" data-i="${i}">Try all ${info.full}</button></div>`;
    if (cls && cls.type === 'osc' && !tr(c.period)) note += `<div class="cat-note warn">It oscillates with period ${cls.period} but has no period set. <button class="small" data-act="setPeriod" data-i="${i}">Set period ${cls.period} and accept any phase</button></div>`;
    if (unstable) note += `<div class="cat-note err">Not stable in this rule, so it can never recover. Remove it or change the rule.</div>`;
    return `<div class="cat${unstable ? ' bad' : ''}" data-i="${i}">
      <canvas width="128" height="128" aria-hidden="true"></canvas>
      <div class="cat-main">
        <div class="cat-top"><span class="name">${esc(catLabel(c, i))}</span>${statusBadge(info)}<span class="dim">${info.w}×${info.h}</span><span class="spacer"></span>
          <button class="small" data-act="dupCat" data-i="${i}">Duplicate</button><button class="small danger" data-act="delCat" data-i="${i}">Remove</button></div>
        <div class="cat-fields">
          <label class="fld"><span class="lbl">Orientations</span><select data-k="cats.${i}.sym" data-re="cats" aria-label="Orientations">${opts}</select></label>
          <label class="fld narrow"><span class="lbl">May be missing for</span><input type="text" data-k="cats.${i}.maxActive" value="${esc(c.maxActive)}" inputmode="numeric" aria-describedby="ma-${i}"><small id="ma-${i}">generations in a row</small></label>
        </div>
        ${note}
        <details class="more"${openSet.has(String(i)) ? ' open' : ''}><summary>More options</summary><div class="body">
          <div class="inline">
            <label class="fld"><span class="lbl">Name</span><input type="text" data-k="cats.${i}.name" value="${esc(c.name)}" placeholder="Catalyst ${i + 1}"></label>
            <label class="fld grow"><span class="lbl">Pattern</span><textarea data-k="cats.${i}.rle" data-rle="1" data-re="cats" rows="2" spellcheck="false">${esc(c.rle)}</textarea></label>
          </div>
          <div class="inline">
            <label class="check" style="align-self:center"><input type="checkbox" data-k="cats.${i}.autoOff" data-re="cats"${c.autoOff ? ' checked' : ''}><span>Centre it automatically <small>offset ${esc(rs.dx)}, ${esc(rs.dy)}</small></span></label>
            ${c.autoOff ? '' : `<label class="fld"><span class="lbl">dx</span><input type="text" data-k="cats.${i}.dx" value="${esc(c.dx)}" inputmode="numeric"></label>
            <label class="fld"><span class="lbl">dy</span><input type="text" data-k="cats.${i}.dy" value="${esc(c.dy)}" inputmode="numeric"></label>`}
            <label class="fld"><span class="lbl">Period <code>period</code></span><input type="text" data-k="cats.${i}.period" value="${esc(c.period)}" inputmode="numeric" placeholder="1"></label>
          </div>
          <div class="checks">
            <label class="check"><input type="checkbox" data-k="cats.${i}.isBlinker"${c.isBlinker ? ' checked' : ''}><span>Accept any phase <code>isblinker</code><small>An oscillator counts as recovered in any of its phases.</small></span></label>
            <label class="check"><input type="checkbox" data-k="cats.${i}.sacrificial" data-re="cats"${c.sacrificial ? ' checked' : ''}><span>Sacrificial <code>sacrificial</code><small>Does not need to recover.</small></span></label>
            <label class="check"><input type="checkbox" data-k="cats.${i}.mustinclude"${c.mustinclude ? ' checked' : ''}><span>Required <code>mustinclude</code><small>Every result uses at least one catalyst marked this way.</small></span></label>
            <label class="check"><input type="checkbox" data-k="cats.${i}.transparent"${c.transparent ? ' checked' : ''}><span>Transparent <code>transparent</code><small>Counts toward "transparent catalysts allowed".</small></span></label>
            <label class="check"><input type="checkbox" data-k="cats.${i}.checkRecovery"${c.checkRecovery ? ' checked' : ''}><span>Check recovery time <code>check-recovery</code><small>It must recover in exactly "may be missing for" generations.</small></span></label>
          </div>
          <h3>Forbidden nearby <span class="opt"><code>forbidden</code>: drop the result if this appears at the catalyst, relative to its top-left</span></h3>
          <div class="sub-rows">${c.forbidden.map((f, j) => subRow(`cats.${i}.forbidden.${j}`, f, 'Forbidden ' + (j + 1), i, j)).join('')}</div>
          <div class="btnrow"><button class="small" data-act="addForbid" data-i="${i}">Add a forbidden pattern</button></div>
          <h3>Cells that must stay on <span class="opt"><code>required</code></span></h3>${subRow(`cats.${i}.req`, c.req, 'Required')}
          <h3>Cells that must stay off <span class="opt"><code>antirequired</code></span></h3>${subRow(`cats.${i}.anti`, c.anti, 'Antirequired')}
          <h3>Cells that must be hit first <span class="opt"><code>locus</code>: or click cells on the thumbnail (red outline)</span></h3>${subRow(`cats.${i}.locus`, c.locus, 'Locus')}
        </div></details>
      </div>
    </div>`;
  }).join('') || '<div class="empty">No catalysts yet. Click a common catalyst above or paste some.</div>';
  $$('#cats .cat').forEach(el => drawCatThumb(el));
  const fix = $('#orientFix'), longs = S.cats.filter(c => +c.maxActive >= 50).length, parts = [];
  if (under) parts.push(`<p><b>${under} catalyst${under > 1 ? 's are' : ' is'} only tried in some orientations.</b> A shape without symmetry has 8 orientations; CatForce only tries the ones its code allows. <button class="small" id="fixAll">Try every orientation for all</button></p>`);
  if (longs) parts.push(`<p><b>${longs} catalyst${longs > 1 ? 's' : ''} may be missing for 50 or more generations,</b> so a catalyst that is destroyed and happens to reform still counts. <button class="small" id="fixMissing">Set all to 20</button></p>`);
  fix.hidden = !parts.length;
  fix.innerHTML = parts.join('');
  renderLib();
}

// Catalyst thumbnail with its locus (cells that must be hit first) in red; click a cell to toggle it.
function thumbGeom(cv, cells) {
  const b = bbox(cells), sc = Math.max(2, Math.min(18, Math.floor((cv.width - 12) / Math.max(b.w, b.h))));
  return {b, sc, ox: Math.floor((cv.width - b.w * sc) / 2), oy: Math.floor((cv.height - b.h * sc) / 2)};
}
function locusCells(c) { const l = c.locus, d = decodeRle(cleanRle(l.rle).rle); return d.cells.map(([x, y]) => [x + (parseInt(l.x, 10) || 0), y + (parseInt(l.y, 10) || 0)]); }
function drawCatThumb(el) {
  const c = S.cats[+el.dataset.i], info = catInfo(cleanRle(c.rle).rle), cv = $('canvas', el);
  drawCells(cv, info.cells);
  if (!info.n) return;
  const g = thumbGeom(cv, info.cells), ctx = cv.getContext('2d');
  ctx.strokeStyle = cssVar('--err'); ctx.lineWidth = 2;
  for (const [x, y] of locusCells(c)) ctx.strokeRect(g.ox + (x - g.b.x0) * g.sc + 1, g.oy + (y - g.b.y0) * g.sc + 1, g.sc - 2, g.sc - 2);
  // required cells (must stay on) in yellow
  const rq = c.req, rd = decodeRle(cleanRle(rq.rle).rle), rx = parseInt(rq.x, 10) || 0, ry = parseInt(rq.y, 10) || 0;
  ctx.strokeStyle = cssVar('--warn');
  for (const [x, y] of rd.cells) ctx.strokeRect(g.ox + (x + rx - g.b.x0) * g.sc + 3, g.oy + (y + ry - g.b.y0) * g.sc + 3, g.sc - 6, g.sc - 6);
}
document.addEventListener('click', e => {
  const cv = e.target.closest && e.target.closest('.cat canvas'); if (!cv) return;
  const el = cv.closest('.cat'), c = S.cats[+el.dataset.i], info = catInfo(cleanRle(c.rle).rle); if (!info.n) return;
  const g = thumbGeom(cv, info.cells), r = cv.getBoundingClientRect();
  const x = Math.floor(((e.clientX - r.left) * cv.width / r.width - g.ox) / g.sc) + g.b.x0, y = Math.floor(((e.clientY - r.top) * cv.height / r.height - g.oy) / g.sc) + g.b.y0;
  if (x < g.b.x0 - 1 || y < g.b.y0 - 1 || x > g.b.x1 + 1 || y > g.b.y1 + 1) return;
  let L = locusCells(c); const k = L.findIndex(p => p[0] === x && p[1] === y);
  if (k >= 0) L.splice(k, 1); else L.push([x, y]);
  if (!L.length) c.locus = {rle: '', x: '0', y: '0'};
  else { const bb = bbox(L); c.locus = {rle: encodeRle(L.map(([a, b2]) => [a - bb.x0, b2 - bb.y0])), x: String(bb.x0), y: String(bb.y0)}; }
  renderCats(); update();
});

// Drag on the map to add a search rectangle (hold Shift to exclude).
(() => {
  const cv = document.getElementById('map'); let st = null;
  const cell = e => { const r = cv.getBoundingClientRect(); return [Math.floor((e.clientX - r.left) / r.width * 64) - 32, Math.floor((e.clientY - r.top) / r.height * 64) - 32]; };
  cv.style.cursor = 'crosshair'; cv.style.touchAction = 'none';
  cv.addEventListener('pointerdown', e => { st = {a: cell(e), ex: e.shiftKey}; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener('pointermove', e => {
    if (!st) return; drawMap(); const b = cell(e), s = cv.width / 64, ctx = cv.getContext('2d');
    ctx.strokeStyle = cssVar(st.ex ? '--err' : '--accent'); ctx.lineWidth = 2;
    ctx.strokeRect((Math.min(st.a[0], b[0]) + 32) * s, (Math.min(st.a[1], b[1]) + 32) * s, (Math.abs(b[0] - st.a[0]) + 1) * s, (Math.abs(b[1] - st.a[1]) + 1) * s);
  });
  cv.addEventListener('pointerup', e => {
    if (!st) return; const b = cell(e), a = st.a;
    S.areas.push({mode: st.ex ? 'exclude' : 'include', x: String(Math.min(a[0], b[0])), y: String(Math.min(a[1], b[1])), w: String(Math.abs(b[0] - a[0]) + 1), h: String(Math.abs(b[1] - a[1]) + 1)});
    st = null; renderAreas(); update();
  });
})();

function renderGoals() {
  const label = {match: 'appears anywhere, in any orientation', andfilter: 'is at this exact spot', orfilter: 'is at this exact spot (or another "one of" goal matches)'};
  $('#goals').innerHTML = S.goals.map((f, i) => `
    <div class="goal" data-i="${i}">
      <label class="fld"><span class="lbl">Goal ${i + 1}: this pattern…</span>
        <select data-k="goals.${i}.type" data-re="goals">${Object.entries(label).map(([v, l]) => `<option value="${v}"${f.type === v ? ' selected' : ''}>${esc(l)}</option>`).join('')}</select></label>
      <label class="fld grow"><span class="lbl">Pattern</span><textarea data-k="goals.${i}.rle" data-rle="1" rows="2" spellcheck="false" placeholder="RLE">${esc(f.rle)}</textarea>
        <span class="chips">${GOAL_CHIPS.map(([n, r]) => `<button class="chip" data-fill="${i}" data-rle="${esc(r)}">${esc(n)}</button>`).join('')}<button class="chip" data-fill="${i}" data-rle="@pat">The active pattern</button></span></label>
      <label class="fld"><span class="lbl">At generation(s)</span><input type="text" data-k="goals.${i}.gen" value="${esc(f.gen)}" placeholder="${esc(defaultGoalGen())}">
        <small>A single number, or a range such as 20-150. A range is only checked once every catalyst has recovered.</small></label>
      ${f.type === 'match' ? '' : `<label class="fld"><span class="lbl">x</span><input type="text" data-k="goals.${i}.dx" value="${esc(f.dx)}" inputmode="numeric"></label>
      <label class="fld"><span class="lbl">y</span><input type="text" data-k="goals.${i}.dy" value="${esc(f.dy)}" inputmode="numeric"></label>`}
      <button class="small danger" data-act="delGoal" data-i="${i}">Remove</button>
    </div>`).join('') || '<div class="empty">No goal: any reaction where the catalysts recover counts.</div>';
}

function renderAreas() {
  $('#areas').innerHTML = S.areas.map((r, i) => `
    <div class="area" data-i="${i}">
      <select data-k="areas.${i}.mode" aria-label="Include or exclude"><option value="include"${r.mode !== 'exclude' ? ' selected' : ''}>Include</option><option value="exclude"${r.mode === 'exclude' ? ' selected' : ''}>Exclude</option></select>
      <input type="text" data-k="areas.${i}.x" value="${esc(r.x)}" inputmode="numeric" placeholder="x" aria-label="x">
      <input type="text" data-k="areas.${i}.y" value="${esc(r.y)}" inputmode="numeric" placeholder="y" aria-label="y">
      <input type="text" data-k="areas.${i}.w" value="${esc(r.w)}" inputmode="numeric" placeholder="width" aria-label="width">
      <input type="text" data-k="areas.${i}.h" value="${esc(r.h)}" inputmode="numeric" placeholder="height" aria-label="height">
      <button class="small danger" data-act="delArea" data-i="${i}">Remove</button>
    </div>`).join('') || '<div class="empty">CatForce default: x −30 to 29, y −30 to 29.</div>';
}

function patCellsAbs() {
  const d = decodeRle(cleanRle(S.patRaw).rle), dx = parseInt(S.patDx, 10) || 0, dy = parseInt(S.patDy, 10) || 0;
  return d.cells.map(([x, y]) => [x + dx, y + dy]);
}
function drawMap() {
  const cv = $('#map'), ctx = cv.getContext('2d'), s = cv.width / 64;
  ctx.clearRect(0, 0, cv.width, cv.height);
  const rects = S.areas.map(r => ({mode: r.mode, v: [r.x, r.y, r.w, r.h].map(t => parseInt(t, 10))})).filter(r => r.v.every(Number.isFinite) && r.v[2] > 0 && r.v[3] > 0);
  const inc = rects.filter(r => r.mode !== 'exclude'), exc = rects.filter(r => r.mode === 'exclude');
  const inR = (r, x, y) => x >= r.v[0] && x < r.v[0] + r.v[2] && y >= r.v[1] && y < r.v[1] + r.v[3];
  ctx.fillStyle = cssVar('--area');
  for (let x = -32; x < 32; x++) for (let y = -32; y < 32; y++) {
    const ok = (inc.length ? inc.some(r => inR(r, x, y)) : (x >= -30 && x < 30 && y >= -30 && y < 30)) && !exc.some(r => inR(r, x, y));
    if (ok) ctx.fillRect((x + 32) * s, (y + 32) * s, s, s);
  }
  ctx.strokeStyle = cssVar('--grid'); ctx.lineWidth = 1;
  for (let i = 0; i <= 64; i += 8) { ctx.beginPath(); ctx.moveTo(i * s + .5, 0); ctx.lineTo(i * s + .5, cv.height); ctx.moveTo(0, i * s + .5); ctx.lineTo(cv.width, i * s + .5); ctx.stroke(); }
  ctx.strokeStyle = cssVar('--axis');
  ctx.beginPath(); ctx.moveTo(32 * s + .5, 0); ctx.lineTo(32 * s + .5, cv.height); ctx.moveTo(0, 32 * s + .5); ctx.lineTo(cv.width, 32 * s + .5); ctx.stroke();
  const pat = patCellsAbs(), key = c => wrap64(c[0]) + ',' + wrap64(c[1]), own = new Set(pat.map(key));
  const sym = CHAIN[S.sym] ? S.sym : 'C1';
  const dot = (x, y) => ctx.fillRect((wrap64(x) + 32) * s, (wrap64(y) + 32) * s, s - 1, s - 1);
  ctx.fillStyle = cssVar('--ghost');
  for (const c of symExpand(pat, sym)) if (!own.has(key(c))) dot(c[0], c[1]);
  ctx.fillStyle = cssVar('--bait');
  for (const b of S.baits) { const d = decodeRle(cleanRle(b.rle).rle), bx = parseInt(b.x, 10) || 0, by = parseInt(b.y, 10) || 0; for (const [x, y] of d.cells) dot(x + bx, y + by); }
  ctx.fillStyle = cssVar('--cell');
  for (const [x, y] of pat) dot(x, y);
  // "Always present" (also-required) cells: yellow outline
  const ad = decodeRle(cleanRle(S.alsoRle).rle), ax = parseInt(S.alsoX, 10) || 0, ay = parseInt(S.alsoY, 10) || 0;
  ctx.strokeStyle = cssVar('--warn'); ctx.lineWidth = 1.5;
  for (const [x, y] of ad.cells) ctx.strokeRect((wrap64(x + ax) + 32) * s + .75, (wrap64(y + ay) + 32) * s + .75, s - 1.5, s - 1.5);
}

function symNote() {
  const box = $('#symMsg'), sym = S.sym;
  if (!sym || sym === 'C1' || !CHAIN[sym]) { box.textContent = ''; return; }
  const pat = patCellsAbs(); if (!pat.length) { box.textContent = ''; return; }
  const k = c => wrap64(c[0]) + ',' + wrap64(c[1]), own = new Set(pat.map(k));
  const all = symExpand(pat, sym), extra = all.filter(c => !own.has(k(c)));
  if (!extra.length) { box.className = 'msg ok'; box.textContent = 'The pattern is already symmetric about this centre.'; return; }
  const touch = extra.some(([x, y]) => { for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) if (own.has(k([x + dx, y + dy]))) return true; return false; });
  box.className = 'msg ' + (touch ? 'warn' : 'tip');
  box.textContent = touch
    ? 'The symmetric copy touches the pattern, so the search starts from the combined ' + all.length + '-cell shape on the map. If that is not what you meant, move the pattern (edit by hand below) or pick another symmetry.'
    : 'The search starts from the pattern plus its symmetric copy (' + all.length + ' cells), shown on the map.';
}

let saveTimer = null;
function update() {
  const r = build(S);
  $('#out').value = r.text;
  const order = {err: 0, warn: 1, tip: 2};
  const list = r.issues.slice().sort((a, b) => order[a.t] - order[b.t]);
  $('#issues').innerHTML = list.length ? list.map(i => `<li class="${i.t}"><span>${esc(i.m)}</span></li>`).join('')
    : '<li class="okline"><span class="ok">Ready to run.</span></li>';
  const rc = checkRule(S.rule), rm = $('#ruleMsg');
  rm.className = 'msg ' + (rc.ok ? 'ok' : 'err');
  rm.textContent = rc.ok ? rc.note[0].toUpperCase() + rc.note.slice(1) + '.' : rc.error[0].toUpperCase() + rc.error.slice(1) + '.';
  $$('[data-rule]').forEach(b => b.classList.toggle('on', b.dataset.rule.toLowerCase() === String(S.rule).replace(/\s+/g, '').toLowerCase()));
  const fn = tr($('#fname').value) || 'search.txt';
  $('#runCheck').textContent = './CatForce ' + fn + ' --check';
  $('#runGo').textContent = './CatForce ' + fn;
  renderPatSummary(); symNote(); drawMap();
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => { try { localStorage.setItem('catforce-gen-v2', JSON.stringify({S, fname: $('#fname').value})); } catch (e) {} }, 400);
}
function renderAll() { syncStatic(); renderObjs(); renderSetups(); renderBaits(); renderCats(); renderGoals(); renderAreas(); update(); }

/* ---------- starting pattern: split a pasted pattern into objects ---------- */
function analysePaste(apply) {
  const msg = $('#pasteMsg'); msg.className = 'msg'; msg.textContent = '';
  const raw = S.paste || '';
  if (!raw.trim()) { split = {objs: [], setups: [], pick: 0}; renderObjs(); renderSetups(); return; }
  const cl = cleanRle(raw);
  if (cl.error) { msg.className = 'msg err'; msg.textContent = 'This does not look like an RLE: ' + cl.error + '.'; return; }
  let ruleNote = '';
  if (cl.rule && apply && cl.rule.replace(/\s/g, '').toLowerCase() !== String(S.rule).replace(/\s/g, '').toLowerCase()) {
    S.rule = cl.rule; $('#f-rule').value = S.rule; catCache = new Map(); renderCats();
    ruleNote = 'Rule set to ' + cl.rule + ' from the pattern header. ';
  }
  const R = ruleOk();
  if (!R.ok) { msg.className = 'msg err'; msg.textContent = 'Fix the rule in step 1 first (' + R.error + ').'; return; }
  const cells = decodeRle(cl.rle).cells;
  if (!cells.length) { msg.className = 'msg err'; msg.textContent = 'The pattern has no live cells.'; return; }
  if (cells.length > 5000) { msg.className = 'msg err'; msg.textContent = 'That pattern is too big to analyse here (' + cells.length + ' cells).'; return; }
  const comps = splitComponents(cells);
  split.objs = comps.map(cc => {
    let info; try { info = classifyObject(new Set(cc.map(p => p.join(','))), R); } catch (e) { info = {type: 'other', label: 'could not be simulated'}; }
    const b = bbox(cc);
    return {cells: cc, info, w: b.w, h: b.h, role: info.type === 'still' || info.type === 'osc' ? 'bait' : 'pattern'};
  });
  if (!split.objs.some(o => o.role === 'pattern')) split.objs[0].role = 'pattern';
  split.pick = 0;
  msg.className = 'msg ok';
  msg.textContent = ruleNote + comps.length + ' object' + (comps.length > 1 ? 's' : '') + ' found. Change a role below if the guess is wrong.';
  recomputeSetups(apply);
}
function recomputeSetups(apply) {
  split.setups = makeSetups(split.objs, Math.max(0, parseInt(S.splitGap, 10) || 0), !!S.splitMerge);
  if (split.pick >= split.setups.length) split.pick = 0;
  if (apply && split.setups.length) applySetup(S, split.setups[split.pick]);
  renderObjs(); renderSetups(); renderBaits(); syncStatic(); update();
}
function applySetup(st, setup) {
  st.patRaw = setup.pat ? setup.pat.rle : '';
  st.patDx = setup.pat ? String(setup.pat.x) : '0';
  st.patDy = setup.pat ? String(setup.pat.y) : '0';
  st.baits = setup.baits.map(b => Object.assign(newBait(), {rle: b.rle, x: String(b.x), y: String(b.y)}));
  st.alsoRle = setup.also ? setup.also.rle : ''; st.alsoX = setup.also ? String(setup.also.x) : '0'; st.alsoY = setup.also ? String(setup.also.y) : '0';
}
function suffixName(name, k) { const m = String(name).match(/^(.*?)(\.[A-Za-z0-9]+)?$/); return m[1] + '_' + k + (m[2] || ''); }
function allSetupFiles() {
  return split.setups.map((st, i) => {
    const S2 = JSON.parse(JSON.stringify(S));
    applySetup(S2, st);
    for (const f of ['output', 'outputList', 'fullReport']) if (tr(S2[f]) || f === 'output') S2[f] = suffixName(tr(S2[f]) || 'results.rle', i + 1);
    const base = tr($('#fname').value) || 'search.txt';
    return {name: suffixName(base, i + 1), text: build(S2).text};
  });
}

/* ---------- catalysts ---------- */
function addCatalyst(rle, name) {
  const info = catInfo(rle);
  if (!info.n) return 'empty';
  const key = canonKey(info.cells);
  if (S.cats.some(c => { const d = decodeRle(cleanRle(c.rle).rle); return d.cells.length && canonKey(d.cells) === key; })) return 'dup';
  const c = newCat(); c.rle = rle; c.name = name || '';
  if (info.cls && info.cls.type === 'osc') { c.period = String(info.cls.period); c.isBlinker = true; }
  S.cats.push(c);
  return 'ok';
}
function pasteCatalysts() {
  const txt = $('#catPaste').value, m = $('#catPasteMsg');
  if (!txt.trim()) { m.className = 'msg err'; m.textContent = 'Paste something first.'; return; }
  let added = 0, dup = 0;
  if (/^\s*cat\s/m.test(txt)) {
    const r = parseIn(txt.split(/\r?\n/).filter(l => /^\s*cat\s/.test(l)).join('\n'));
    for (const c of r.st.cats) { S.cats.push(c); added++; }
  } else {
    const cl = cleanRle(txt);
    if (cl.error) { m.className = 'msg err'; m.textContent = cl.error; return; }
    const comps = splitComponents(decodeRle(cl.rle).cells);
    for (const g of clusterObjects(comps.map(bbox), 1)) {
      const cells = [].concat(...g.map(k => comps[k])), b = bbox(cells);
      const r = addCatalyst(encodeRle(cells.map(([x, y]) => [x - b.x0, y - b.y0])));
      if (r === 'ok') added++; else if (r === 'dup') dup++;
    }
  }
  m.className = 'msg ' + (added ? 'ok' : 'warn');
  m.textContent = 'Added ' + added + ' catalyst' + (added === 1 ? '' : 's') + (dup ? '; ' + dup + ' repeat' + (dup > 1 ? 's' : '') + ' skipped' : '') + '.';
  if (added) $('#catPaste').value = '';
  renderCats(); update();
}

/* ---------- loading ---------- */
function loadState(st, opts = {}) {
  S = Object.assign(defaultState(), st);
  split = {objs: [], setups: [], pick: 0};
  catCache = new Map();
  $('#exampleNote').hidden = !opts.example;
  renderAll();
  if (S.paste) analysePaste(false);
}
const EXAMPLE = (() => {
  const st = defaultState();
  Object.assign(st, {
    rule: 'B3/S23', paste: 'x = 3, y = 4, rule = B3/S23\no$3o$obo$2bo!', patRaw: 'o$3o$obo$2bo!', patDx: '-1', patDy: '-2',
    maxGen: '100', numCat: '2', stable: '15', output: 'herschel.rle', maxTime: '2h',
    goals: [{type: 'match', gen: '20-100', rle: 'o$3o$obo$2bo!', dx: '0', dy: '0'}]
  });
  for (const n of ['Block', 'Eater 1', 'Beehive', 'Tub', 'Boat', 'Snake']) {
    const c = newCat(); c.name = n; c.rle = LIB.find(l => l[0] === n)[1]; st.cats.push(c);
  }
  return st;
})();

/* ---------- events ---------- */
let pasteTimer = null;
document.addEventListener('input', e => {
  const el = e.target;
  if (el.id === 'f-includes') { S.includes = el.value.split('\n'); update(); return; }
  if (el.id === 'fname') { update(); return; }
  const k = el.dataset && el.dataset.k; if (!k) return;
  setPath(S, k, el.type === 'checkbox' ? el.checked : el.value);
  if (k === 'paste') { clearTimeout(pasteTimer); pasteTimer = setTimeout(() => analysePaste(true), 350); return; }
  if (k === 'rule') { catCache = new Map(); renderCats(); }
  if (k === 'splitGap' || k === 'splitMerge') { recomputeSetups(true); return; }
  if (el.type === 'checkbox' || el.tagName === 'SELECT') {
    if (el.dataset.re === 'cats') renderCats();
    if (el.dataset.re === 'goals') renderGoals();
  }
  update();
});
document.addEventListener('change', e => {
  const el = e.target;
  if (el.classList.contains('role')) { split.objs[+el.closest('.obj').dataset.k].role = el.value; recomputeSetups(true); return; }
  if (el.name === 'setupPick') { split.pick = +el.value; applySetup(S, split.setups[split.pick]); renderSetups(); renderBaits(); syncStatic(); update(); return; }
  if (el.dataset && el.dataset.rle && el.dataset.k) {      // tidy a pasted Golly RLE once editing finishes
    const c = cleanRle(el.value);
    if (c.rle && !c.error && c.rle !== el.value) { el.value = c.rle; setPath(S, el.dataset.k, c.rle); }
    if (el.dataset.re === 'cats') renderCats();
    update();
  }
});
document.addEventListener('click', async e => {
  const b = e.target.closest('button'); if (!b) return;
  const i = b.dataset.i !== undefined ? +b.dataset.i : -1;
  if (b.dataset.rule) { S.rule = b.dataset.rule; $('#f-rule').value = S.rule; catCache = new Map(); renderCats(); update(); return; }
  if (b.dataset.lib !== undefined) { const [n, r] = LIB[+b.dataset.lib]; addCatalyst(r, n); renderCats(); update(); return; }
  if (b.dataset.quad) { const q = b.dataset.quad.split(' '); S.areas.push({mode: 'exclude', x: q[0], y: q[1], w: q[2], h: q[3]}); renderAreas(); update(); return; }
  if (b.dataset.fill !== undefined) {
    e.preventDefault();
    const gi = +b.dataset.fill;
    if (b.dataset.rle === '@pat') {
      const all = symExpand(patCellsAbs(), CHAIN[S.sym] ? S.sym : 'C1');
      if (!all.length) return;
      const bb = bbox(all);
      S.goals[gi].rle = encodeRle(all.map(([x, y]) => [x - bb.x0, y - bb.y0])); S.goals[gi].dx = String(bb.x0); S.goals[gi].dy = String(bb.y0);
    } else S.goals[gi].rle = b.dataset.rle;
    renderGoals(); update(); return;
  }
  switch (b.dataset.act || b.id) {
    case 'delCat': S.cats.splice(i, 1); renderCats(); break;
    case 'dupCat': S.cats.splice(i + 1, 0, JSON.parse(JSON.stringify(S.cats[i]))); renderCats(); break;
    case 'autoSym': S.cats[i].sym = 'auto'; renderCats(); break;
    case 'setPeriod': { const info = catInfo(cleanRle(S.cats[i].rle).rle); S.cats[i].period = String(info.cls.period); S.cats[i].isBlinker = true; renderCats(); break; }
    case 'fixAll': S.cats.forEach(c => { c.sym = 'auto'; }); renderCats(); break;
    case 'fixMissing': S.cats.forEach(c => { if (+c.maxActive >= 50) c.maxActive = '20'; }); renderCats(); break;
    case 'addForbid': S.cats[i].forbidden.push({rle: '', x: '0', y: '0'}); renderCats(); break;
    case 'delForbid': S.cats[i].forbidden.splice(+b.dataset.j, 1); renderCats(); break;
    case 'delBait': S.baits.splice(i, 1); renderBaits(); break;
    case 'addBait': S.baits.push(newBait()); renderBaits(); break;
    case 'delGoal': S.goals.splice(i, 1); renderGoals(); break;
    case 'addMatch': S.goals.push(newGoal('match')); renderGoals(); break;
    case 'addAnd': S.goals.push(newGoal('andfilter')); renderGoals(); break;
    case 'addOr': S.goals.push(newGoal('orfilter')); renderGoals(); break;
    case 'delArea': S.areas.splice(i, 1); renderAreas(); break;
    case 'addArea': S.areas.push(newArea('include')); renderAreas(); break;
    case 'addExclude': S.areas.push(newArea('exclude')); renderAreas(); break;
    case 'catPasteGo': pasteCatalysts(); return;
    case 'btnOpen': $('#openPanel').hidden = false; $('#importBox').focus(); return;
    case 'importClose': $('#openPanel').hidden = true; return;
    case 'importAll': {
      const r = parseIn($('#importBox').value);
      if (!r.st.patRaw && !r.st.cats.length) { $('#importMsg').className = 'msg err'; $('#importMsg').textContent = 'No pat or cat lines found in that text.'; return; }
      loadState(r.st);
      $('#importMsg').className = 'msg ' + (r.ignored.length ? 'warn' : 'ok');
      $('#importMsg').textContent = 'Loaded ' + r.st.cats.length + ' catalysts and ' + r.st.goals.length + ' goals.' + (r.ignored.length ? ' Skipped: ' + r.ignored.slice(0, 3).join('; ') + (r.ignored.length > 3 ? '…' : '') : '');
      return;
    }
    case 'importCats': {
      const r = parseIn($('#importBox').value.split(/\r?\n/).filter(l => /^\s*cat\s/.test(l)).join('\n'));
      S.cats.push(...r.st.cats); renderCats(); update();
      $('#importMsg').className = 'msg ok'; $('#importMsg').textContent = 'Added ' + r.st.cats.length + ' catalysts.';
      return;
    }
    case 'btnExample': loadState(JSON.parse(JSON.stringify(EXAMPLE)), {example: true}); return;
    case 'btnReset': $('#resetConfirm').hidden = false; return;
    case 'resetNo': $('#resetConfirm').hidden = true; return;
    case 'resetYes': $('#resetConfirm').hidden = true; loadState(defaultState()); return;
    case 'copy': copyText($('#out').value, $('#copyMsg'), $('#out')); return;
    case 'save': saveFiles([{name: tr($('#fname').value) || 'search.txt', text: $('#out').value}], $('#copyMsg')); return;
    case 'copyAll': {
      const files = allSetupFiles();
      copyText(files.map(f => '# ===== ' + f.name + ' =====\n' + f.text).join('\n'), $('#allMsg'), $('#out'));
      return;
    }
    case 'saveAll': saveFiles(allSetupFiles(), $('#allMsg')); return;
    default: return;
  }
  update();
});
$('#importFile').addEventListener('change', e => {
  const f = e.target.files && e.target.files[0]; if (!f) return;
  const rd = new FileReader();
  rd.onload = () => { $('#importBox').value = String(rd.result); $('#importMsg').className = 'msg ok'; $('#importMsg').textContent = 'Read ' + f.name + '. Choose what to load.'; };
  rd.readAsText(f);
});

async function copyText(text, msgEl, fallbackEl) {
  msgEl.className = 'msg ok';
  try { await navigator.clipboard.writeText(text); msgEl.textContent = 'Copied.'; return; } catch (e) {}
  fallbackEl.value = text; fallbackEl.focus(); fallbackEl.select();
  msgEl.className = 'msg warn'; msgEl.textContent = 'Selected: press Ctrl+C or ⌘C to copy.';
}
let dl = null;
(async () => {
  try { if (window.claude && claude.use) dl = await claude.use('downloads'); } catch (e) { dl = null; }
  if (dl) { $('#save').hidden = false; const sa = $('#saveAll'); if (sa) sa.hidden = false; }
})();
async function saveFiles(files, msgEl) {
  if (!dl) return;
  try {
    for (const f of files) await dl.save({filename: f.name.replace(/[^\w.\-]+/g, '_'), data: f.text});
    msgEl.className = 'msg ok'; msgEl.textContent = files.length > 1 ? 'Saved ' + files.length + ' files.' : 'Saved ' + files[0].name + '.';
  } catch (err) {
    if (!err || err.code !== 'declined') { msgEl.className = 'msg warn'; msgEl.textContent = 'Could not save; use Copy instead.'; }
  }
}

/* ---------- start ---------- */
initStatic();
(() => {
  let saved = null;
  try { saved = JSON.parse(localStorage.getItem('catforce-gen-v2') || 'null'); } catch (e) { saved = null; }
  if (saved && saved.S) { if (saved.fname) $('#fname').value = saved.fname; loadState(saved.S); }
  else loadState(JSON.parse(JSON.stringify(EXAMPLE)), {example: true});
})();
if (window.matchMedia) matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { renderObjs(); renderCats(); update(); });

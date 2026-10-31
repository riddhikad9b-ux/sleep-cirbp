# # Sleep = the brain's nightly "bowhead whale mode"?
# **Hypothesis.** When we sleep, the brain cools by 1–2 °C. Cooling switches on **CIRBP**, the cold-inducible protein that gives the 200-year-old bowhead whale its better DNA repair (Firsanov et al., *Nature* 2025). Sleep loss keeps the brain warm, so CIRBP falls (it is the No. 1 down gene after sleep loss in mouse hippocampus, Gaine et al. 2021) and DNA double-strand breaks pile up in neurons (Bellesi 2016; Zada 2019).
# 
# **This notebook tests 4 predictions on public data**
# 
# | Test | Prediction if hypothesis is right | Data |
# |---|---|---|
# | T1 | *Cirbp* (and its partner *Rbm3*) go **down** after sleep loss, in more than one dataset, and come back after recovery sleep | GSE166831, GSE33302 (hippocampus); GSE113754, GSE140345 (cortex) |
# | T2 | Inside each dataset, animals with more *Cirbp* also show more DNA-repair gene activity | same |
# | T3 | People with too little sleep show lower *CIRBP* in blood | GSE39445 (human, 6 h vs 10 h sleep) |
# | T4 | Exercise (which warms the body) does **not** raise *Cirbp*; it may lower it | GSE330623 (rat), GSE9565 (mouse) |
# 
# Run: **Runtime → Run all**. At the end copy the **SEND THIS TO CLAUDE** box.
# A failed test prints FAILED with the reason and the rest still runs.

!pip -q install GEOparse 2>/dev/null | grep -v -i dependency || true
import os, re, requests, traceback, warnings
import numpy as np, pandas as pd
import GEOparse
from scipy import stats
warnings.filterwarnings('ignore'); os.makedirs('data', exist_ok=True)
REPORT = []
def note(s): print(s); REPORT.append(str(s))

NHEJ = ['XRCC5','XRCC6','PRKDC','LIG4','XRCC4','NHEJ1','DCLRE1C','PAXX','POLL','POLM']
HR   = ['RAD51','BRCA1','BRCA2','RAD51C','RAD51B','MRE11','MRE11A','RAD50','NBN','RBBP8','PALB2','EXO1','BLM']
DDR  = ['ATM','H2AX','H2AFX','MDC1','TP53BP1','PARP1','RNF8','RNF168']
COLD = ['CIRBP','RBM3']
POS  = ['HOMER1','ARC','FOS','EGR1','NR4A1']   # known to RISE after sleep loss: checks the labels are right
ORIGCOL = {}                                  # dataset -> {GSM id: original column name in the count file}
import time, glob
def get_text(url, tries=4):
    for i in range(tries):
        try:
            r = requests.get(url, timeout=90); r.raise_for_status(); return r.text
        except Exception as e:
            print('  retry', i + 1, url, e); time.sleep(5 * (i + 1))
    raise IOError(f'could not reach {url}')
def get_geo(gid, **kw):
    for i in range(4):
        try: return GEOparse.get_GEO(geo=gid, destdir='data', silent=True, **kw)
        except Exception as e:
            print('  GEO retry', i + 1, gid, e)
            for f in glob.glob(f'data/{gid}*'): os.remove(f)    # remove the half-downloaded file
            time.sleep(10 * (i + 1))
    raise IOError(f'{gid}: download kept failing (NCBI busy); run this cell again later')
REPAIR = NHEJ + HR + DDR

# ## Helpers: load a GEO dataset, find groups, turn probes into genes

def norm(s): return s.astype(str).str.lower().str.replace(r'[^a-z0-9]+', ' ', regex=True).str.strip()

def varying_text(meta):
    cols = [c for c in meta.columns if c == 'title' or c.startswith('characteristics')
            or c.startswith('source_name') or c.startswith('description')]
    keep = [c for c in cols if meta[c].astype(str).nunique() > 1] or ['title']
    return meta[keep].apply(norm).agg(' | '.join, axis=1)

def label(meta, case_re, ctrl_re, drop_re=None, gid=''):
    t = varying_text(meta); g = pd.Series('other', index=t.index)
    c, k = t.str.contains(case_re), t.str.contains(ctrl_re)
    g[c & ~k] = 'case'; g[k & ~c] = 'control'
    g[c & k & t.str.contains(r'\bnon |\bnsd\b|\bno |\bun')] = 'control'   # 'non sleep deprived' = control
    if drop_re: g[t.str.contains(drop_re)] = 'other'
    vc = g.value_counts().to_dict()
    note(f'{gid}: case={vc.get("case",0)} control={vc.get("control",0)} other={vc.get("other",0)}')
    s = pd.DataFrame({'g': g, 't': t.str.replace(r'\b(rep|replicate|sample)\s*\d*\b|\b\d+\b', '', regex=True)
                      .str.replace(r'\s+', ' ', regex=True).str.slice(0, 70)})
    for (gg, tt), n in s.value_counts().sort_index().head(12).items(): note(f'    {gg:8s} x{n:<3d} {tt}')
    return g

SYM_RE = r'^[A-Za-z][A-Za-z0-9\-\.]{1,19}$'
def first_symbol(v, col):
    if not isinstance(v, str) or v.strip() in ('', '---', 'nan', 'None'): return None
    parts = [p.strip() for p in v.split(' /// ')[0].split(' // ')]
    if col.lower() == 'gene_assignment': return parts[1] if len(parts) > 1 else None
    return parts[0] or None

def probes_to_genes(expr, gpl):
    gpl = gpl.copy(); gpl['ID'] = gpl['ID'].astype(str)
    names = ['gene symbol','gene_symbol','symbol','gene symbols','gene_assignment','gene list','genesymbol','gene_name','name']
    sym = None
    for col in [c for c in gpl.columns if c.lower() in names]:
        s = gpl.set_index('ID')[col].apply(lambda v: first_symbol(v, col)).dropna()
        s = s[s.str.match(SYM_RE)]
        if len(s) > 0.3 * len(gpl): sym = s; break
    if sym is None: raise ValueError('no gene-name column: ' + ', '.join(gpl.columns[:12]))
    e = expr.copy(); e.index = e.index.astype(str); e = e.apply(pd.to_numeric, errors='coerce')
    e = e[e.index.isin(sym.index)]; e['gene'] = sym.reindex(e.index).str.upper().values
    e['m'] = e.drop(columns='gene').mean(1)
    g = e.dropna(subset=['m']).sort_values('m', ascending=False).drop_duplicates('gene').set_index('gene').drop(columns='m')
    g = g.dropna(thresh=max(2, int(0.8 * g.shape[1]))); g = g.apply(lambda r: r.fillna(r.mean()), axis=1)
    if g.max().max() > 50: g = np.log2(g.clip(lower=1))
    return g

def load_array(gid):
    gse = get_geo(gid)
    meta = gse.phenotype_data
    pid = meta['platform_id'].mode()[0] if 'platform_id' in meta else list(gse.gpls)[0]
    meta = meta[meta.get('platform_id', pid) == pid] if 'platform_id' in meta else meta
    expr = gse.pivot_samples('VALUE')[meta.index]
    return meta, probes_to_genes(expr, gse.gpls[pid].table)

def key(s): return re.sub(r'[^a-z0-9]', '', str(s).lower())

def map_ids(c, species, gid):
    # rows = Ensembl or NCBI gene IDs -> gene symbols; rescues CIRBP/RBM3 if the bulk lookup missed them
    os.system('pip -q install mygene'); import mygene; mg = mygene.MyGeneInfo()
    idx = c.index.astype(str).str.split('.').str[0]
    if idx.str.upper().str.startswith('ENS').mean() > 0.5: scope = 'ensembl.gene'
    elif idx.str.fullmatch(r'\d+').mean() > 0.5: scope = 'entrezgene'
    else:
        c = c.copy(); c.index = idx.str.upper(); return c.groupby(level=0).sum()
    c = c.copy(); c.index = idx
    q = mg.querymany(list(idx.unique()), scopes=scope, fields='symbol', species=species, as_dataframe=True, verbose=False)
    if 'symbol' not in q.columns:
        note(f'{gid}: !! ID lookup found nothing; first IDs: {list(idx[:5])}')
        sym = pd.Series(dtype=str)
    else:
        sym = q['symbol'].dropna().astype(str).str.upper(); sym = sym[~sym.index.duplicated()]
    for g in COLD + ['PARP1']:
        if g in set(sym.values): continue
        hits = mg.query(f'symbol:{g.capitalize()}', species=species, fields='ensembl.gene,entrezgene').get('hits', [])
        note(f'{gid}: {g} not in bulk lookup; direct search found IDs: ' +
             ', '.join(str(h.get('entrezgene', '')) + '/' + str(h.get('ensembl', '')) for h in hits[:3])[:300])
        for h in hits:
            ids = [str(h.get('entrezgene', ''))]
            e = h.get('ensembl', [])
            for x in (e if isinstance(e, list) else [e]):
                gg = x.get('gene') if isinstance(x, dict) else None
                ids += gg if isinstance(gg, list) else [gg]
            for i in ids:
                if i and i in c.index: sym[i] = g; note(f'{gid}: rescued {g} from ID {i}')
    c = c[c.index.isin(sym.index)]; c.index = sym.reindex(c.index).values
    note(f'{gid}: {c.index.nunique()} genes after mapping {scope} -> symbol ({species}); CIRBP present: {"CIRBP" in set(c.index)}')
    return c.groupby(level=0).sum()

def load_counts(gid, species='mouse'):
    # authors' count file -> log2 CPM, gene symbols
    meta = get_geo(gid, how='brief').phenotype_data
    try:   # GEO's own count table first (NCBI gene IDs)
        page = get_text(f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={gid}')
        links = re.findall(r'(/geo/download/\?type=rnaseq_counts[^"\']*raw_counts[^"\']*)', page)
        if links:
            c = pd.read_csv('https://www.ncbi.nlm.nih.gov' + links[0].replace('&amp;', '&'), sep='\t', index_col=0, compression='gzip')
            os.system('pip -q install mygene'); import mygene
            c = map_ids(c, species, gid)
            note(f'{gid}: GEO count table, {c.shape[0]} genes x {c.shape[1]} samples')
            logcpm = np.log2(c.div(c.sum(0), axis=1) * 1e6 + 1)
            return meta, logcpm[(c >= 10).sum(1) >= 3]
    except Exception as e:
        print(gid, 'GEO count table failed:', e)
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gid[:-3]}nnn/{gid}/suppl/'
    KNOWN = {'GSE166831': 'GSE166831_counts_all.csv.gz', 'GSE330623': 'GSE330623_processed_data_raw.csv.gz',
             'GSE140345': 'GSE140345_RNA_level_all_genes.txt.gz', 'GSE113754': 'GSE113754_RNASeq_S3_WT_SD5_HC5_Counts.txt.gz'}
    try: files = [f for f in re.findall(r'href="([^"]+)"', get_text(base)) if gid in f and not f.endswith('/')]
    except Exception as e: files = []
    if not any(re.search(r'count|raw|processed|rna_level', f, re.I) for f in files) and gid in KNOWN:
        files = [KNOWN[gid]]; note(f'{gid}: file list unavailable, using known file name')
    pick = [f for f in files if re.search(r'count|raw|processed|rna|genes|expression', f, re.I)
            and not re.search(r'atac|peak|chip', f, re.I) and re.search(r'\.(txt|tsv|csv)(\.gz)?$', f, re.I)]
    note(f'{gid}: author file used: {pick[:1]}')
    if not pick: raise ValueError(f'no count file in {files}')
    f = pick[0]; raw = pd.read_csv(base + f, sep=',' if '.csv' in f.lower() else '\t',
                                   compression='gzip' if f.endswith('.gz') else None)
    txt = [c for c in raw.columns if not pd.api.types.is_numeric_dtype(raw[c])]
    sc = next((c for c in txt if re.search(r'symbol|gene_?name|^name$|^gene$', str(c), re.I)), None)
    c = raw.set_index(sc or (txt[0] if txt else raw.columns[0])).select_dtypes('number')
    c = c[[x for x in c.columns if not re.search(r'length|^start|^end|^chr|gc', str(x), re.I)]]
    if not set(c.columns) & set(meta.index):
        tk = {key(t): g for g, t in meta['title'].items()}; m = {}
        for col in c.columns:
            k = key(col); hit = tk.get(k) or next((g for t, g in tk.items() if k and (k in t or t in k)), None)
            if hit and hit not in m.values(): m[col] = hit
        how = 'matched to sample titles'
        if len(m) < c.shape[1] and c.shape[1] == len(meta):
            m = dict(zip(c.columns, meta.index)); how = 'FILE ORDER (titles did not match; check labels)'
        c = c[list(m)].rename(columns=m)
        ORIGCOL[gid] = {v: k for k, v in m.items()}
        note(f'{gid}: columns {how}; e.g. ' + ', '.join(f'{k} -> {v}' for k, v in list(m.items())[:3]))
    c.index = c.index.astype(str).str.upper()
    if not sc and c.index.str.contains(r'[|_:]').mean() > 0.5 and c.index.str.contains('ENS').mean() > 0.5:
        c.index = c.index.str.extract(r'(ENS[A-Z]*G\d+)')[0].fillna(c.index.to_series().reset_index(drop=True)).values
        note(f'{gid}: pulled Ensembl IDs out of combined row names')
    note(f'{gid}: first row names: {list(c.index[:4])}')
    c = c.groupby(level=0).sum()
    if not sc: c = map_ids(c, species, gid)
    logcpm = np.log2(c.div(c.sum(0), axis=1) * 1e6 + 1)
    return meta, logcpm[(c >= 10).sum(1) >= 3]

# ## Statistics used
# - **Effect:** log2 fold change (case − control) and **Hedges' g** (standardised difference, bias-corrected).
# - **p-value:** Welch t-test. **Rank:** where *Cirbp* sits among all genes (1 = most changed).
# - **T2 co-expression:** within each group separately (so the sleep-loss effect itself can't create the link), Spearman correlation of *Cirbp* with a **repair score** (mean z-score of NHEJ + HR genes), and a **null** from 1000 random gene sets of the same size.

def hedges_g(a, b):
    na, nb = len(a), len(b); s = np.sqrt(((na-1)*np.var(a, ddof=1) + (nb-1)*np.var(b, ddof=1)) / (na+nb-2))
    d = (np.mean(a) - np.mean(b)) / s if s > 0 else np.nan
    return d * (1 - 3 / (4*(na+nb) - 9)), np.sqrt((na+nb)/(na*nb) + d**2/(2*(na+nb)))

def compare(expr, grp, gid, genes=COLD + POS + ['PARP1', 'XRCC5', 'XRCC6', 'PRKDC', 'RAD51']):
    a = expr[[s for s in grp[grp == 'case'].index if s in expr.columns]]
    b = expr[[s for s in grp[grp == 'control'].index if s in expr.columns]]
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    res = pd.DataFrame({'log2FC': a.mean(1) - b.mean(1), 't': t, 'p': p}, index=expr.index).dropna()
    res['rank'] = res['p'].rank()
    rows = []
    for gname in genes:
        if gname in res.index:
            gg, se = hedges_g(a.loc[gname].values, b.loc[gname].values)
            r = res.loc[gname]
            rows.append(dict(gene=gname, log2FC=round(r.log2FC, 3), g=round(gg, 2), se=round(se, 2),
                             p=f'{r.p:.2g}', rank=f'{int(r["rank"])}/{len(res)}'))
    out = pd.DataFrame(rows, columns=['gene', 'log2FC', 'g', 'se', 'p', 'rank'])
    if out.empty:
        note(f'{gid}: !! none of {genes[:3]} found. First gene names: {list(expr.index[:8])}')
    else:
        note(f'{gid} (n case={a.shape[1]}, control={b.shape[1]}):\n' + out.to_string(index=False))
        pc = out[out.gene.isin(POS)]
        if len(pc):
            up = (pc.log2FC > 0).sum()
            note(f'    label check: {up}/{len(pc)} known sleep-loss genes go UP '
                 f'({"labels look right" if up >= max(3, len(pc) - 1) else "!! labels may be wrong or effect weak"})')
    return res, out

def repair_coexp(expr, grp, gid, n_null=1000, seed=0):
    if 'CIRBP' not in expr.index: note(f'{gid}: CIRBP not measured'); return
    rep = [g for g in NHEJ + HR if g in expr.index]
    rng = np.random.default_rng(seed); pool = [g for g in expr.index if g not in rep + COLD]
    for lev in ['control', 'case']:
        s = [x for x in grp[grp == lev].index if x in expr.columns]
        if len(s) < 6: note(f'{gid} {lev}: only {len(s)} samples, skip'); continue
        z = expr[s].sub(expr[s].mean(1), axis=0).div(expr[s].std(1).replace(0, np.nan), axis=0)
        rho = stats.spearmanr(z.loc['CIRBP'], z.loc[rep].mean()).correlation
        null = [stats.spearmanr(z.loc['CIRBP'], z.loc[list(rng.choice(pool, len(rep), replace=False))].mean()).correlation
                for _ in range(n_null)]
        null = np.array(null); p = (np.sum(null >= rho) + 1) / (n_null + 1)
        note(f'{gid} {lev} (n={len(s)}): CIRBP vs repair score rho={rho:.2f}; '
             f'random gene sets median rho={np.nanmedian(null):.2f}; p(higher than random)={p:.3f}; {len(rep)} repair genes')

# ## T1 + T2. Mouse hippocampus after sleep loss

T1, RECOVERY = {}, {}
try:
    meta, ex = load_counts('GSE166831')
    grp = label(meta, r'^sd(?:\b|\d)|\| sd\b', r'^nsd(?:\b|\d)|\| nsd\b', gid='GSE166831')
    if (grp == 'case').sum() < 2:   # titles are 'NSD'/'SD'
        title = norm(meta['title']); grp = pd.Series('other', index=meta.index)
        grp[title.str.match(r'^sd(?:\b|\d)')] = 'case'; grp[title.str.match(r'^nsd(?:\b|\d)')] = 'control'
        note(f'GSE166831 title rule: {grp.value_counts().to_dict()}')
    T1['GSE166831'] = compare(ex, grp, 'T1 GSE166831 hippocampus RNA-seq, 5 h sleep loss')
    repair_coexp(ex, grp, 'T2 GSE166831')
except Exception as e:
    note(f'GSE166831 FAILED: {type(e).__name__}: {e}')

try:
    meta, ex = load_array('GSE33302')
    grp = label(meta, r'sleep ?depriv|\bsd\b|deprived', r'control|\bnsd\b|non ?sleep|undisturbed|home ?cage', r'recover', gid='GSE33302')
    T1['GSE33302'] = compare(ex, grp, 'T1 GSE33302 hippocampus microarray')
    repair_coexp(ex, grp, 'T2 GSE33302')
except Exception as e:
    note(f'GSE33302 FAILED: {type(e).__name__}: {e}')

# ## T1b. Cortex after sleep loss, and after recovery sleep (other labs)
# - **GSE113754:** frontal cortex, 5 h sleep loss vs home-cage controls, 5 vs 5.
# - **GSE140345:** cortex, SD 3 h and 6 h vs time-matched home-cage controls, plus 6 h recovery sleep after SD (RS6) vs its control (HC12).
# 
# Groups are matched **by time of day** (SD3 vs HC3, SD6 vs HC6, RS6 vs HC12), because *Cirbp* also follows the daily clock.
# **Extra prediction:** after recovery sleep, *Cirbp* should come back to normal (RS6 ≈ HC12).

def timepoint_groups(meta, gid=None):
    t = norm(meta['title']) + ' | ' + varying_text(meta)
    if gid in ORIGCOL:   # the count file's own column names say which group each column is
        oc = norm(pd.Series({k: ORIGCOL[gid].get(k, '') for k in meta.index}))
        a = oc.str.extract(r'\b(?:wt|s3)?(sd|hc|rs|nsd)\s*(\d{1,2})\b'); b = t.str.extract(r'\b(?:wt|s3)?(sd|hc|rs|nsd)\s*(\d{1,2})\b')
        both = a[0].notna() & b[0].notna()
        bad = ((a[0] + a[1]) != (b[0] + b[1]))[both].sum()
        if bad: note(f'{gid}: !! {bad} columns: file name and GEO title disagree; using the FILE column names')
        t = oc + ' | ' + t      # file column name first, so its group tag wins
    tag = t.str.extract(r'\b(?:wt|s3)?(sd|hc|rs|nsd)\s*(\d{1,2})\b')
    tp = (tag[0] + tag[1]).where(tag[0].notna())
    alt = t.str.extract(r'\bt\s*(\d{1,2})\s+(nsd|sd|rs|hc)\b')     # e.g. 'rna t6 sd x': time first, then group
    tp = tp.fillna((alt[1].replace({'nsd': 'hc'}) + alt[0]).where(alt[0].notna()))
    return tp, t

for gid, pairs in [('GSE113754', [('sd5', 'hc5')]),
                   ('GSE140345', 'auto')]:
    try:
        meta, ex = load_counts(gid)
        tp, t = timepoint_groups(meta, gid)
        mut = t.str.split(' | ', regex=False).str[0].str.contains(r'^s3|shank') if gid in ORIGCOL else t.str.contains('shank3')
        if mut.any():     # GSE113754 has Shank3 mutants: keep wild-type only (judged from the file column name)
            tp = tp.where(~mut)
            note(f'{gid}: Shank3 mutant samples left out')
        note(f'{gid} time groups found: {tp.value_counts().to_dict()}')
        if tp.notna().sum() == 0:
            note('    sample text examples: ' + ' || '.join(t.drop_duplicates().head(4)))
        if pairs == 'auto':   # every time point that has both SD (or RS) and a control
            have = set(tp.dropna())
            pairs = [(g, 'hc' + g[2:]) for g in sorted(have) if g[:2] in ('sd', 'rs') and 'hc' + g[2:] in have]
            note(f'{gid}: comparisons found: {pairs}')
        for case, ctrl in pairs:
            grp = pd.Series('other', index=meta.index); grp[tp == case] = 'case'; grp[tp == ctrl] = 'control'
            if (grp == 'case').sum() < 2 or (grp == 'control').sum() < 2:
                if gid == 'GSE113754':   # fall back to plain SD vs HC labels
                    grp = label(meta, r'\bsd\b|sleep ?depriv', r'\bhc\b|home ?cage|control', r'recover|\brs\b', gid=gid)
                else:
                    note(f'{gid} {case} vs {ctrl}: not enough samples'); continue
            res = compare(ex, grp, f'T1b {gid} cortex {case.upper()} vs {ctrl.upper()}')
            if case.startswith('rs'): RECOVERY[gid] = res
            elif case.startswith('sd'):
                T1[gid] = res     # one comparison per dataset in the meta-analysis: the last (longest) SD time
    except Exception as e:
        note(f'{gid} FAILED: {type(e).__name__}: {e}')

# ### T1 combined: random-effects meta-analysis of *Cirbp*
# DerSimonian–Laird pooling of Hedges' g across datasets that loaded.

def meta_g(gs, ses):
    gs, ses = np.array(gs, float), np.array(ses, float); w = 1 / ses**2
    fe = np.sum(w*gs) / np.sum(w); q = np.sum(w*(gs-fe)**2); dfq = len(gs) - 1
    tau2 = max(0, (q - dfq) / (np.sum(w) - np.sum(w**2)/np.sum(w))) if dfq > 0 else 0
    wr = 1 / (ses**2 + tau2); pooled = np.sum(wr*gs) / np.sum(wr); se = np.sqrt(1/np.sum(wr))
    p = 2 * stats.norm.sf(abs(pooled/se)); return pooled, se, p, tau2
for gene in COLD:
    rows = [o[1].set_index('gene').loc[gene] for o in T1.values()
            if isinstance(o[1], pd.DataFrame) and 'gene' in o[1].columns and gene in set(o[1]['gene'])]
    note(f'T1 META {gene}: {len(rows)} dataset(s) with this gene')
    if len(rows) >= 2:
        pg, se, p, tau2 = meta_g([r.g for r in rows], [r.se for r in rows])
        note(f'T1 META {gene}: pooled g={pg:.2f} (95% CI {pg-1.96*se:.2f} to {pg+1.96*se:.2f}), p={p:.2g}, '
             f'datasets={len(rows)}, tau2={tau2:.2f}')

# ## T3. Humans: 1 week of 6 h sleep vs 10 h sleep (blood)
# GSE39445 (Möller-Levet et al., *PNAS* 2013). Many blood samples per person across a day, so we average each person's samples per condition first (a paired test, person by person).

try:
    meta, ex = load_array('GSE39445')
    t = varying_text(meta)
    note('GSE39445 example sample text: ' + ' || '.join(t.drop_duplicates().head(4)))
    cond = pd.Series('other', index=t.index)
    cond[t.str.contains(r'sleep restriction|restrict|insufficient')] = 'short'
    cond[t.str.contains(r'sleep extension|extension|sufficient') & ~t.str.contains(r'restrict|insufficient')] = 'long'
    subj = t.str.extract(r'subject\s+([a-z]{0,3}\d+)')[0]
    note(f'GSE39445 conditions: {cond.value_counts().to_dict()}, people found: {subj.nunique()}')
    df = pd.DataFrame({'cond': cond, 'subj': subj}).dropna()
    df = df[df.cond != 'other']
    for gene in COLD + ['PARP1', 'XRCC5', 'XRCC6']:
        if gene not in ex.index: continue
        v = ex.loc[gene, df.index]
        per = pd.DataFrame({'v': v, 'cond': df.cond, 'subj': df.subj}).groupby(['subj', 'cond'])['v'].mean().unstack()
        per = per.dropna()
        if len(per) >= 5:
            d = per['short'] - per['long']; w = stats.wilcoxon(d)
            note(f'T3 {gene}: short-minus-long sleep = {d.mean():+.3f} log2 (people={len(per)}), '
                 f'{(d < 0).sum()}/{len(d)} people lower after short sleep, Wilcoxon p={w.pvalue:.2g}')
        else:
            note(f'T3 {gene}: could not pair people (found {len(per)}); check labels above')
except Exception as e:
    note(f'GSE39445 FAILED: {type(e).__name__}: {e}')

# ## T4. Exercise: does it raise or lower *Cirbp*?
# Exercise heats the body. If *Cirbp* tracks temperature, regular running should **not** raise it. This tells us whether "exercise reverses sleep loss" can hold for this gene.

try:
    meta, ex = load_counts('GSE330623', species='rat')
    title = norm(meta['title']); grp = pd.Series('other', index=meta.index)
    grp[title.str.contains('exercise control')] = 'case'; grp[title.str.contains('sedentary control')] = 'control'
    note(f'GSE330623 groups: {grp.value_counts().to_dict()}')
    compare(ex, grp, 'T4 GSE330623 rat hippocampus, exercise vs sedentary')
except Exception as e:
    note(f'GSE330623 FAILED: {type(e).__name__}: {e}')
try:
    meta, ex = load_array('GSE9565')
    grp = label(meta, r'exercis|runn|wheel|\brun\b', r'sedentary|control|\bsed\b|locked', gid='GSE9565')
    ch = [c for c in meta.columns if re.search(r'source_name|label_ch|characteristics', c)]
    note('GSE9565 channel fields: ' + ' || '.join(f'{c}={meta[c].iloc[0]}' for c in ch[:6]))
    if (grp == 'case').sum() >= 2 and (grp == 'control').sum() >= 2:
        compare(ex, grp, 'T4 GSE9565 mouse hippocampus, running vs sedentary')
    else:   # two-colour arrays: each sample is already a log ratio of the two channels
        note(f'GSE9565: Cy5 = exercise, Cy3 = sedentary, so a positive log ratio = higher with exercise. '
             f'CIRBP on array: {"CIRBP" in ex.index}; genes on array: {len(ex)}')
        for gene in COLD + POS + ['PARP1']:
            if gene in ex.index:
                v = ex.loc[gene].values; tt = stats.ttest_1samp(v, 0)
                note(f'T4 GSE9565 {gene}: mean log ratio {v.mean():+.3f} over {len(v)} arrays, p={tt.pvalue:.2g} '
                     '(two-colour; sign depends on which channel is exercise, see channel fields above)')
except Exception as e:
    note(f'GSE9565 FAILED: {type(e).__name__}: {e}')

# ## How to read the result
# - **Supports the idea:** T1 pooled g clearly below 0; T2 rho above the random sets in controls; T3 most people lower after short sleep.
# - **Weakens it:** T2 no link to repair genes (then *Cirbp* may only be a thermometer, not part of the repair story); T3 no change in humans.
# - **T4** is a side question: if exercise lowers *Cirbp*, then for this gene exercise does **not** undo sleep loss; you need the sleep itself.

print('=' * 70); print('SEND THIS TO CLAUDE'); print('=' * 70); print('\n'.join(REPORT)); print('=' * 70)

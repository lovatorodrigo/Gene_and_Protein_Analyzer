# app_streamlit.py
# UI (Streamlit) — execução via SUBPROCESS (shim)
# + Auditoria BLAST, presets, histórico de runs, favoritos, highlight multicor com tooltip,
#   densidade de colchetes, bundle zip, KPIs por fonte, links rápidos, mapeador proteína→códon,
#   CATÁLOGO DE PADRÕES (por categorias) + RNAi (siRNA/shRNA) + geração de miRNA seed (7mer/8mer) + GC%.

import os
import sys
import io
import re
import json
import shutil
import zipfile
import traceback
import subprocess
from pathlib import Path
from datetime import datetime

import streamlit as st

# ==== Config inicial ====
st.set_page_config(page_title="Effatha — Functional Regions", layout="wide")

# Optional deps
try:
    import pandas as pd
except Exception:
    pd = None

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

# -------------------------------
# CSS — wrap visual, multicor e tooltip (alto contraste + fix p/ spans curtos)
# -------------------------------
st.markdown(
    """
    <style>
    .stCode > div > pre, .stCode pre, pre, code {
        white-space: pre-wrap !important;
        word-break: break-word !important;
        overflow-wrap: anywhere !important;
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace !important;
        font-size: 12.5px !important;
        line-height: 1.4 !important;
    }
    .mono-pre {
        white-space: pre-wrap !important;
        word-break: break-word !important;
        overflow-wrap: anywhere !important;
        font-family: ui-monospace, monospace !important;
        font-size: 12.5px !important;
        border: 1px solid #e6e6e6;
        border-radius: 6px;
        padding: 10px;
        background: #ffffff;
        overflow: visible;              /* evita clip de tooltip */
        position: relative;
    }

    /* Badges/KPIs */
    .hl { background: #fff2a8; padding: 0 1px; border-radius: 2px; }
    .badge { display:inline-block; padding:1px 6px; border-radius: 10px; font-size: 11px; margin-right: 6px; }
    .badge.full { background:#e7f1ff; }
    .badge.uniprot { background:#e6ffe7; }
    .badge.genpept { background:#fff1e6; }
    .badge.lig { background:#ffe6f2; }
    .kpibox { background:#f7f7f9; padding:10px; border-radius:8px; border:1px solid #eee; margin-bottom:8px;}

    /* Highlights (alto contraste) */
    .hl-tip {
        position: relative;
        cursor: help;
        border-radius: 3px;
        padding: 0 2px;
        box-shadow: inset 0 -2px 0 rgba(0,0,0,0.18);
        font-weight: 600;
    }
    .hl-aa    { background:#FFC107; color:#2b1e00; border-bottom:2px solid #d39e00; }
    .hl-nt    { background:#29B6F6; color:#002231; border-bottom:2px solid #0288D1; }
    .hl-rnai  { background:#66BB6A; color:#06230a; border-bottom:2px solid #2E7D32; }
    .hl-prom  { background:#BA68C8; color:#230028; border-bottom:2px solid #8E24AA; }
    .hl-motif { background:#FF7043; color:#311000; border-bottom:2px solid #E64A19; }
    .hl-custom{ background:#FFD54F; color:#332400; border-bottom:2px solid #FFA000; }

    /* Tooltip centralizado, com largura mínima (fix p/ matches curtos) */
    .hl-tip:hover::after {
        content: attr(data-title);
        position: absolute;
        left: 50%;
        top: 100%;
        transform: translate(-50%, 8px);
        background: #222;
        color: #fff;
        padding: 8px 10px;
        border-radius: 4px;
        font-size: 12px;
        line-height: 1.25;
        z-index: 99999;
        box-shadow: 0 4px 12px rgba(0,0,0,0.35);
        /* >>> as três linhas abaixo resolvem o “efeito coluna” */
        width: max-content;                                  /* usa largura do conteúdo */
        min-width: 220px;                                     /* nunca fica estreito */
        max-width: min(560px, calc(100vw - 56px));            /* não passa da viewport */
        white-space: pre-wrap;                                /* quebra em espaços/\\n */
        text-align: left;
        pointer-events: none;
    }
    .hl-tip:hover::before {
        content: "";
        position: absolute;
        left: 50%;
        top: 100%;
        transform: translate(-50%, 0);
        border: 6px solid transparent;
        border-top-color: #222;
        z-index: 99998;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# -------------------------------
# Utils básicos
# -------------------------------
def make_run_dir(base_dir="runs"):
    d = (Path(base_dir) / datetime.now().strftime("%Y%m%d_%H%M%S")).resolve()
    d.mkdir(parents=True, exist_ok=True)
    return d

def pick_card(cards, tag, start, end):
    for c in cards or []:
        if c.get("tag")==tag and c.get("start")==start and c.get("end")==end:
            return c
    return None

def deep_update(d, u):
    for k, v in (u or {}).items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            deep_update(d[k], v)
        else:
            d[k] = v

def write_shim(run_dir, module_path, cfg_overrides):
    cfg_path = run_dir / "cfg_overrides.json"
    cfg_path.write_text(json.dumps(cfg_overrides, indent=2, ensure_ascii=False), encoding="utf-8")

    shim_code = f"""# AUTO-GERADO pelo app Streamlit (não editar)
import json, sys, importlib.util, pathlib, traceback
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

def deep_update(d, u):
    for k, v in (u or {{}}).items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            deep_update(d[k], v)
        else:
            d[k] = v

module_path = r\"\"\"{Path(module_path).resolve()}\"\"\" 
spec = importlib.util.spec_from_file_location("pipeline_mod", module_path)
if spec is None or spec.loader is None:
    raise RuntimeError("Falha ao carregar spec para o pipeline: " + module_path)

mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod  # importante p/ dataclasses/typing
try:
    spec.loader.exec_module(mod)  # executa o main.py
except Exception:
    traceback.print_exc()
    raise

with open("cfg_overrides.json", "r", encoding="utf-8") as fh:
    overrides = json.load(fh)
deep_update(mod.CONFIG, overrides)

src = (mod.CONFIG["input"]["source"] or "").lower()
if src == "uniprot":
    mod.run_from_uniprot(mod.CONFIG["input"]["uniprot_acc"])
elif src == "pdb":
    mod.run_from_pdb(mod.CONFIG["input"]["pdb_id"])
elif src == "ncbi_protein":
    mod.run_from_ncbi_protein(mod.CONFIG["input"]["ncbi_protein_acc"])
elif src == "ncbi_gene":
    mod.run_from_ncbi_gene(mod.CONFIG["input"]["gene"])
else:
    raise ValueError("CONFIG['input']['source'] deve ser 'pdb' | 'uniprot' | 'ncbi_protein' | 'ncbi_gene'.")
"""
    shim_path = run_dir / "shim_run.py"
    shim_path.write_text(shim_code, encoding="utf-8")
    return shim_path

def run_subprocess(py_exe, script_path, cwd):
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(
        [py_exe, script_path.name],
        cwd=str(cwd),
        text=True,
        capture_output=True,
        env=env
    )
    return proc.returncode, proc.stdout, proc.stderr

# ---- History & Presets ----
def list_runs(base_dir="runs", max_items=20):
    base = Path(base_dir)
    if not base.exists():
        return []
    items = [p for p in base.iterdir() if p.is_dir()]
    items.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return items[:max_items]

def load_context(run_dir: Path):
    ctx_path = run_dir / "context_summary.json"
    if ctx_path.exists():
        try:
            return json.loads(ctx_path.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None

def build_quick_links(protein_meta: dict):
    acc = (protein_meta or {}).get("accession") or ""
    links = []
    if acc:
        # UniProt ACC
        if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9](-\d+)?$", acc) or re.match(r"^[A-NR-Z]\d{5}(-\d+)?$", acc):
            links.append(("UniProt", f"https://www.uniprot.org/uniprotkb/{acc}/entry"))
        # RefSeq Protein ACC
        if re.match(r"^[NXYWP]P_\d+(?:\.\d+)?$", acc, re.I):
            links.append(("NCBI Protein", f"https://www.ncbi.nlm.nih.gov/protein/{acc}"))
    return links

# -------------------------------
# Busca/realce: suporte multicor + tooltip + catálogo
# -------------------------------
def html_escape(s: str) -> str:
    return (s
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;"))

def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq.upper() if c in "GC")
    return 100.0 * gc / len(seq)

def revcomp_rna(s: str) -> str:
    t = s.upper().replace("T","U")
    comp = {"A":"U","U":"A","G":"C","C":"G"}
    return "".join(comp.get(c,"N") for c in t[::-1])

def revcomp_dna(s: str) -> str:
    t = s.upper().replace("U","T")
    comp = {"A":"T","T":"A","G":"C","C":"G"}
    return "".join(comp.get(c,"N") for c in t[::-1])

# Esquema de cores por categoria
CATEGORY_CLASS = {
    "AA": "hl-aa",
    "NT": "hl-nt",
    "RNAi": "hl-rnai",
    "PROM": "hl-prom",
    "MOTIF": "hl-motif",
    "CUSTOM": "hl-custom",
}

# Catálogo de padrões (por categoria)
# Cada entrada: dict(
#   key, label, regex, flags, scope ('AA'|'DNA'|'mRNA'|'NT'), category_tag (p/ cor),
#   tooltip (str|callable(match_str)->str), priority (int pequeno = maior prioridade)
# )
PATTERN_CATALOG = {
    "Nucleotídeos — Microssatélites / repetições": [
        {"key":"NT_AT_GE10_DNA", "label":"[AT]{10,} (DNA)", "regex":r"[AT]{10,}", "flags":re.I, "scope":"DNA", "category_tag":"NT", "tooltip":"Microssatélite AT≥10 (DNA)", "priority":30},
        {"key":"NT_AU_GE10_mRNA","label":"[AU]{10,} (mRNA)","regex":r"[AU]{10,}", "flags":re.I, "scope":"mRNA","category_tag":"NT","tooltip":"Microssatélite AU≥10 (mRNA)","priority":30},
        {"key":"NT_CG_RUN_DNA",  "label":"CGCGCG (≥3 rep.)", "regex":r"(?:CG){3,}", "flags":0,   "scope":"DNA", "category_tag":"NT", "tooltip":"Dinucleotídeo CG repetido (DNA)", "priority":40},
    ],
    "Nucleotídeos — Start/Stop (DNA; use U para mRNA)": [
        {"key":"NT_ATG_START", "label":"ATG (start)","regex":r"ATG","flags":0,"scope":"DNA","category_tag":"PROM","tooltip":"Códon inicial (DNA)","priority":50},
        {"key":"NT_STOP_TAA",  "label":"TAA (stop)","regex":r"TAA","flags":0,"scope":"DNA","category_tag":"PROM","tooltip":"Códon stop (DNA)","priority":60},
        {"key":"NT_STOP_TAG",  "label":"TAG (stop)","regex":r"TAG","flags":0,"scope":"DNA","category_tag":"PROM","tooltip":"Códon stop (DNA)","priority":60},
        {"key":"NT_STOP_TGA",  "label":"TGA (stop)","regex":r"TGA","flags":0,"scope":"DNA","category_tag":"PROM","tooltip":"Códon stop (DNA)","priority":60},
    ],
    "Nucleotídeos — Promotores e boxes": [
        {"key":"NT_TATA_BOX","label":"TATA box (~-25)","regex":r"TATA[AT]A[AT]","flags":re.I,"scope":"DNA","category_tag":"PROM","tooltip":"TATA box (DNA)","priority":70},
        {"key":"NT_CAAT_BOX","label":"CAAT box","regex":r"CAAT","flags":re.I,"scope":"DNA","category_tag":"PROM","tooltip":"CAAT box (DNA)","priority":70},
        {"key":"NT_GC_BOX","label":"GC box (SP1)","regex":r"GGGCGG","flags":re.I,"scope":"DNA","category_tag":"PROM","tooltip":"GC box (SP1) (DNA)","priority":70},
    ],
    "RNAi (siRNA/shRNA)": [
        # siRNA candidatos (DNA/mRNA) com GC% no núcleo de 19nt
        {
            "key":"RNAI_SIRNA_DNA_AA_N19_TT",
            "label":"siRNA DNA: AA[ACGT]{19}TT",
            "regex":r"AA[ACGT]{19}TT",
            "flags":re.I,
            "scope":"DNA",
            "category_tag":"RNAi",
            "tooltip":lambda s: f"RNAi candidato (DNA) — núcleo 19nt GC={gc_percent(s[2:-2]):.1f}%",
            "priority":10
        },
        {
            "key":"RNAI_SIRNA_mRNA_AA_N19_UU",
            "label":"siRNA mRNA: AA[ACGU]{19}UU",
            "regex":r"AA[ACGU]{19}UU",
            "flags":re.I,
            "scope":"mRNA",
            "category_tag":"RNAi",
            "tooltip":lambda s: f"RNAi candidato (mRNA) — núcleo 19nt GC={gc_percent(s[2:-2]):.1f}%",
            "priority":10
        },
        # shRNA loop e terminadores
        {"key":"RNAI_SHRNA_LOOP_DNA","label":"Loop shRNA (DNA) TTCAAGAGA","regex":r"TTCAAGAGA","flags":0,"scope":"DNA","category_tag":"RNAi","tooltip":"Loop shRNA (DNA) TTCAAGAGA","priority":20},
        {"key":"RNAI_SHRNA_LOOP_mRNA","label":"Loop shRNA (mRNA) UUCAAGAGA","regex":r"UUCAAGAGA","flags":0,"scope":"mRNA","category_tag":"RNAi","tooltip":"Loop shRNA (mRNA) UUCAAGAGA","priority":20},
        {"key":"RNAI_POLIII_TERM_DNA","label":"Pol III terminator T{4,}","regex":r"T{4,}","flags":0,"scope":"DNA","category_tag":"RNAi","tooltip":"Pol III terminator — evitar em hairpin (DNA)","priority":25},
        {"key":"RNAI_POLIII_TERM_mRNA","label":"Pol III terminator U{4,}","regex":r"U{4,}","flags":0,"scope":"mRNA","category_tag":"RNAi","tooltip":"Pol III terminator — evitar em hairpin (mRNA)","priority":25},
        {"key":"RNAI_ARE_DNA","label":"ARE candidato ATTTTA (DNA)","regex":r"ATTTTA","flags":0,"scope":"DNA","category_tag":"RNAi","tooltip":"AU-rich element (DNA proxy)","priority":35},
        {"key":"RNAI_ARE_mRNA","label":"ARE candidato AUUUA (mRNA)","regex":r"AUUUA","flags":0,"scope":"mRNA","category_tag":"RNAi","tooltip":"AU-rich element (mRNA)","priority":35},
    ],
    "Aminoácidos — Sinalização / localização": [
        {"key":"AA_N_GLYCO", "label":"N-glicosilação N[^P][ST][^P]", "regex":r"N[^P][ST][^P]", "flags":0, "scope":"AA", "category_tag":"MOTIF", "tooltip":"Motivo N-glicosilação", "priority":30},
        {"key":"AA_NLS_BASIC","label":"NLS básico (K/R cluster)","regex":r"[KR]{3,}","flags":0,"scope":"AA","category_tag":"MOTIF","tooltip":"Sinal nuclear básico (heurístico)","priority":40},
        {"key":"AA_PKA","label":"PKA [RK][RK]x[ST]","regex":r"[RK]{2}.{1}[ST]","flags":0,"scope":"AA","category_tag":"MOTIF","tooltip":"Motivo fosforilação PKA (consenso simpl.)","priority":40},
        {"key":"AA_SH3","label":"SH3 PxxP","regex":r"P..P","flags":0,"scope":"AA","category_tag":"MOTIF","tooltip":"Motivo prólina PxxP (SH3-like)","priority":50},
        {"key":"AA_ZIPPER","label":"Leu zipper Lx6Lx6Lx6L","regex":r"L.{6}L.{6}L.{6}L","flags":0,"scope":"AA","category_tag":"MOTIF","tooltip":"Leucine zipper (aprox.)","priority":60},
    ],
    "Aminoácidos — Metais / enzimas": [
        {"key":"AA_C2H2","label":"Zn-finger C2H2 (aprox.)","regex":r"C.{2,4}C.{12}H.{3,5}H","flags":0,"scope":"AA","category_tag":"MOTIF","tooltip":"C2H2 Zn-finger (aprox.)","priority":60},
        {"key":"AA_HIS_CLUSTER","label":"His cluster H{3,}","regex":r"H{3,}","flags":0,"scope":"AA","category_tag":"MOTIF","tooltip":"Cluster de histidinas (ligantes/ácidos)","priority":60},
    ],
}

def pattern_specs_for_scope(all_specs, scope_key):
    """Filtra specs por escopo ('AA', 'DNA', 'mRNA')."""
    out = []
    for sp in all_specs:
        sc = sp.get("scope")
        if sc == scope_key or (scope_key in ("DNA","mRNA") and sc=="NT"):
            out.append(sp)
    return out

def compile_specs(specs):
    """Compila regex e organiza prioridade."""
    out = []
    for sp in specs:
        try:
            rx = re.compile(sp["regex"], sp.get("flags", 0))
            out.append({**sp, "_rx": rx})
        except Exception:
            # ignora padrão inválido
            continue
    out.sort(key=lambda d: int(d.get("priority", 100)))
    return out

def apply_multi_highlight(text: str, specs, color_override=None):
    """
    Realce multicor com tooltip; evita sobrepor matches (primeiro que encaixa ganha).
    specs: lista de dicts compilados com campos:
       - _rx (regex), label, category_tag, tooltip (str | callable(match_str)->str)
    """
    if not text or not specs:
        return f'<div class="mono-pre">{html_escape(text or "")}</div>'

    # Vamos marcar índices ocupados para evitar overlap
    taken = [False] * len(text)
    intervals = []  # (start, end, class, tooltip, label)

    for sp in specs:
        rx = sp.get("_rx")
        if not rx:
            continue
        cat = sp.get("category_tag","CUSTOM")
        css_class = CATEGORY_CLASS.get(cat, CATEGORY_CLASS["CUSTOM"])
        lab = sp.get("label","padrão")
        tooltip = sp.get("tooltip")
        for m in rx.finditer(text):
            s, e = m.span()
            if s == e:
                continue
            # skip se já houver cobertura
            if any(taken[i] for i in range(s, e)):
                continue
            # calcula tooltip
            tip = lab
            if callable(tooltip):
                try:
                    tip = tooltip(m.group(0))
                except Exception:
                    tip = lab
            elif isinstance(tooltip, str) and tooltip:
                tip = tooltip
            # ocupa
            for i in range(s, e):
                taken[i] = True
            intervals.append((s, e, css_class, tip, lab))

    # Monta HTML
    intervals.sort(key=lambda x: x[0])
    out = []
    cur = 0
    for s, e, klass, tip, lab in intervals:
        # trecho "limpo"
        if cur < s:
            out.append(html_escape(text[cur:s]))
        # trecho realçado
        frag = html_escape(text[s:e])
        out.append(f'<span class="hl-tip {klass}" data-title="{html_attr(tip)}">{frag}</span>')
        cur = e
    if cur < len(text):
        out.append(html_escape(text[cur:]))

    return f'<div class="mono-pre">{"".join(out)}</div>'

# ---- Densidade de colchetes Effatha ----
def effatha_variant_flags(seq: str):
    flags = []
    i = 0
    L = len(seq or "")
    while i < L:
        ch = seq[i]
        if ch == "[":
            j = seq.find("]", i+1)
            if j == -1:
                flags.append(0)
                i += 1
            else:
                flags.append(1)
                i = j + 1
        elif ch == "*" and i == L - 1:
            i += 1
        else:
            flags.append(0)
            i += 1
    return flags

def bracket_density(seq: str):
    flags = effatha_variant_flags(seq or "")
    if not flags:
        return 0.0, []
    d = sum(flags) / len(flags)
    return d, flags

def moving_avg(vals, k=15):
    if not vals:
        return []
    out = []
    n = len(vals)
    for i in range(n):
        a = max(0, i - k//2)
        b = min(n, i + k//2 + 1)
        out.append(sum(vals[a:b]) / (b - a))
    return out

# ---- Favorites ----
def fav_key(tag, start, end):
    return f"{tag}:{start}-{end}"

def is_fav(tag, start, end):
    fset = st.session_state.get("favorites", set())
    return fav_key(tag, start, end) in fset

def toggle_fav(tag, start, end, v):
    fset = set(st.session_state.get("favorites", set()))
    key = fav_key(tag, start, end)
    if v:
        fset.add(key)
    else:
        fset.discard(key)
    st.session_state["favorites"] = fset

# ---- BLAST audit helpers ----
def guess_accession_from_sbjct(s: str):
    if not s:
        return None, None
    m = re.search(r"\|([A-Z0-9]{1,6}\d{1,5}(?:-\d+)?)\|", s)
    if m:
        acc = m.group(1)
        return "uniprot", acc
    m = re.search(r"(?:^|[|\\s])(NP|XP|YP|WP)_\\d+(?:\\.\\d+)?", s)
    if m:
        acc = m.group(0).strip(" |")
        return "refseq_protein", acc
    if re.match(r"^[A-NR-Z]\\d{5}(-\\d+)?$", s) or re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9](-\\d+)?$", s):
        return "uniprot", s
    return None, None

def make_ext_link(kind, acc):
    if not kind or not acc:
        return ""
    if kind == "uniprot":
        return f"https://www.uniprot.org/uniprotkb/{acc}/entry"
    if kind == "refseq_protein":
        return f"https://www.ncbi.nlm.nih.gov/protein/{acc}"
    return ""

def html_attr(s) -> str:
    """Escapa string para uso em atributos HTML (ex.: data-title="...")."""
    if s is None:
        return ""
    s = str(s)
    return (s.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace('"', "&quot;")
             .replace("'", "&#39;")
             .replace("\n", "\\n"))

# -------------------------------
# Estado inicial
# -------------------------------
st.title("Gene & Protein Analyzer")

if "favorites" not in st.session_state:
    st.session_state["favorites"] = set()

if "expand_all" not in st.session_state:
    st.session_state["expand_all"] = False

if "last_loaded_run" not in st.session_state:
    st.session_state["last_loaded_run"] = None

# -------------------------------
# Sidebar — Configurações / Presets / Histórico / Avançado
# -------------------------------
with st.sidebar:
    st.header("⚙️ Configurações")

    # Histórico de runs
    st.subheader("Histórico de execuções")
    base_runs_dir = st.text_input("Pasta de runs", value="runs", key="base_runs_dir")

    def _list_runs(base_dir="runs", max_items=30):
        base = Path(base_dir)
        if not base.exists():
            return []
        items = [p for p in base.iterdir() if p.is_dir()]
        items.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return items[:max_items]

    runs = _list_runs(base_runs_dir, max_items=30)
    run_labels = [f"{p.name}  ({datetime.fromtimestamp(p.stat().st_mtime).strftime('%Y-%m-%d %H:%M:%S')})" for p in runs]
    idx_load = st.selectbox("Carregar run anterior", options=["(nenhum)"] + run_labels, index=0)
    if idx_load != "(nenhum)":
        sel = runs[run_labels.index(idx_load)]
        ctx = load_context(sel)
        if ctx:
            st.session_state["last_run_dir"] = str(sel)
            st.session_state["last_context"] = ctx
            st.session_state["last_paths"] = {
                "report": str(next(sel.glob("*.txt"), sel / "report.txt")),
                "regions_csv": str(next(sel.glob("regions_*.csv"), sel / "regions.csv")),
                "blast_csv": str(next(sel.glob("variants_blast_*.csv"), "")),
                "prompt_md": str(next(sel.glob("prompt_effatha_*.md"), "")),
                "prompt_txt": str(next(sel.glob("prompt_effatha_*.txt"), "")),
            }
            st.session_state["last_enable_blast"] = True
            st.session_state["last_loaded_run"] = sel.name
            st.success(f"Run {sel.name} carregado do histórico.")
        else:
            st.warning("context_summary.json não encontrado nesse run.")

    st.markdown("---")

    if "source_choice" not in st.session_state:
        st.session_state["source_choice"] = "UniProt"

    pipeline_path = st.text_input("Caminho do pipeline (main.py)",
                                  value=str(Path("main.py").resolve()),
                                  key="pipeline_path")

    c1, c2 = st.columns([3,1])
    with c1:
        st.radio("Fonte de entrada",
                 options=["UniProt", "PDB", "NCBI Protein", "NCBI Gene"],
                 horizontal=True, key="source_choice")
    with c2:
        if st.button("Limpar entradas & resultados"):
            for k in ("uniprot_acc","pdb_id","ncbi_protein_acc","gene_id","gene_symbol"):
                if k in st.session_state:
                    st.session_state[k] = ""
            for k in ("last_run_dir","last_context","last_paths","last_enable_blast","last_loaded_run"):
                if k in st.session_state:
                    st.session_state[k] = None
            st.session_state["favorites"] = set()
            st.rerun()

    st.text_input("UniProt ACC", value="P02458", key="uniprot_acc")
    st.text_input("PDB ID", value="5xnl", key="pdb_id")
    st.text_input("NCBI Protein (RefSeq)", value="NP_000537.3", key="ncbi_protein_acc")

    st.markdown("**NCBI Gene**")
    colg1, colg2, colg3 = st.columns([1,1,1])
    with colg1:
        st.radio("ID", options=["GeneID","Símbolo"], horizontal=True, index=0, key="gene_mode")
    with colg2:
        st.text_input("GeneID", value="7157", key="gene_id")
        st.text_input("Símbolo", value="TP53", key="gene_symbol")
    with colg3:
        st.number_input("TaxID (p/ Símbolo)", min_value=1, value=9606, key="gene_taxid")
        st.selectbox("Isoforma", options=["longest","mane"], index=0, key="isoform_policy")

    st.markdown("---")
    st.subheader("Camadas")
    st.checkbox("UniProt VARIANT", value=True, key="enable_uniprot_variant")
    st.checkbox("Proteins Variation API", value=True, key="enable_variation_api")
    st.checkbox("BLAST", value=False, key="enable_blast")

    st.markdown("---")
    st.subheader("Filtros de variação (linha FILTRADA)")
    st.checkbox("Ativar filtros", value=True, key="enable_filters")
    col1, col2 = st.columns(2)
    with col1:
        st.number_input("Frequência alélica mínima (AF)", min_value=0.0, max_value=1.0,
                        step=0.001, value=0.01, key="min_af")
        st.selectbox("Evidência (ECO)", options=["prefer_manual","manual_only","any"], index=0, key="evidence_mode")
        st.selectbox("Evidência ausente", options=["include","exclude"], index=0, key="missing_evid")
    with col2:
        st.multiselect("Clínica permitida",
                       options=["benign","likely_benign","association","protective"],
                       default=["benign","likely_benign"], key="allow_clin")
        st.checkbox("Permitir 'conflicting' se presente", value=False, key="allow_if_conflicting")
        st.checkbox("Permitir 'uncertain' se sem rótulo", value=False, key="allow_uncertain")

    st.caption("Predições: manter variantes preditas benignas/toleradas.")
    colp = st.columns(3)
    with colp[0]:
        st.multiselect("SIFT allow",
                       options=["tolerated","tolerated_low_confidence","deleterious"],
                       default=["tolerated","tolerated_low_confidence"], key="sift_allow")
    with colp[1]:
        st.multiselect("PolyPhen allow",
                       options=["benign","possibly_benign","possibly_damaging","damaging"],
                       default=["benign","possibly_benign"], key="polyphen_allow")
    with colp[2]:
        st.checkbox("Requerer SIFT e PolyPhen OK", value=False, key="require_both")

    st.checkbox("Incluir BLAST na linha FILTRADA", value=False,
                disabled=not st.session_state.get("enable_blast", False),
                key="include_blast_in_filtered")

    st.markdown("---")
    st.subheader("Avançado (BLAST/Regiões)")
    with st.expander("Opções avançadas", expanded=False):
        st.number_input("BLASTp: min identidade (proporção)", min_value=0.5, max_value=1.0, value=0.97, step=0.01, key="blast_min_id")
        st.number_input("BLASTp: min cobertura (proporção)", min_value=0.5, max_value=1.0, value=0.95, step=0.01, key="blast_min_cov")
        st.selectbox("BLASTp: db", options=["swissprot","nr"], index=0, key="blast_db")
        st.number_input("BLASTp: hitlist_size", min_value=1, max_value=200, value=25, step=1, key="blast_hitlist_size")
        st.checkbox("Restringir por espécie (entrez_query por taxid)", value=True, key="blast_same_species_only")

        st.markdown("---")
        st.number_input("Flanco pontual (AA) ±", min_value=0, max_value=50, value=5, step=1, key="point_flank")
        st.number_input("Comprimento mínimo de região", min_value=1, max_value=200, value=6, step=1, key="default_min_len")
        st.checkbox("Mesclar overlaps (mesmo tipo)", value=True, key="merge_overlaps")
        st.checkbox("Incluir FULL 1..L", value=True, key="add_full")
        st.checkbox("Incluir features NCBI (GenPept) quando não mapear UniProt", value=True, key="include_ncbi_features")

    st.markdown("---")
    st.subheader("Presets")
    presets_dir = Path(base_runs_dir) / "presets"
    presets_dir.mkdir(parents=True, exist_ok=True)
    preset_name = st.text_input("Nome do preset", value="meu_preset.json", key="preset_name")
    if st.button("💾 Salvar preset"):
        preset = {
            "source_choice": st.session_state.get("source_choice"),
            "uniprot_acc": st.session_state.get("uniprot_acc"),
            "pdb_id": st.session_state.get("pdb_id"),
            "ncbi_protein_acc": st.session_state.get("ncbi_protein_acc"),
            "gene_mode": st.session_state.get("gene_mode"),
            "gene_id": st.session_state.get("gene_id"),
            "gene_symbol": st.session_state.get("gene_symbol"),
            "gene_taxid": st.session_state.get("gene_taxid"),
            "isoform_policy": st.session_state.get("isoform_policy"),
            "enable_uniprot_variant": st.session_state.get("enable_uniprot_variant"),
            "enable_variation_api": st.session_state.get("enable_variation_api"),
            "enable_blast": st.session_state.get("enable_blast"),
            "enable_filters": st.session_state.get("enable_filters"),
            "min_af": st.session_state.get("min_af"),
            "evidence_mode": st.session_state.get("evidence_mode"),
            "missing_evid": st.session_state.get("missing_evid"),
            "allow_clin": st.session_state.get("allow_clin"),
            "allow_if_conflicting": st.session_state.get("allow_if_conflicting"),
            "allow_uncertain": st.session_state.get("allow_uncertain"),
            "sift_allow": st.session_state.get("sift_allow"),
            "polyphen_allow": st.session_state.get("polyphen_allow"),
            "require_both": st.session_state.get("require_both"),
            "include_blast_in_filtered": st.session_state.get("include_blast_in_filtered"),
            "point_flank": st.session_state.get("point_flank"),
            "default_min_len": st.session_state.get("default_min_len"),
            "merge_overlaps": st.session_state.get("merge_overlaps"),
            "add_full": st.session_state.get("add_full"),
            "include_ncbi_features": st.session_state.get("include_ncbi_features"),
            "blast_min_id": st.session_state.get("blast_min_id"),
            "blast_min_cov": st.session_state.get("blast_min_cov"),
            "blast_db": st.session_state.get("blast_db"),
            "blast_hitlist_size": st.session_state.get("blast_hitlist_size"),
            "blast_same_species_only": st.session_state.get("blast_same_species_only"),
        }
        (presets_dir / preset_name).write_text(json.dumps(preset, indent=2, ensure_ascii=False), encoding="utf-8")
        st.success(f"Preset salvo em {presets_dir / preset_name}")

    preset_files = [p.name for p in presets_dir.glob("*.json")]
    sel_preset = st.selectbox("Carregar preset", options=["(nenhum)"] + preset_files, index=0)
    if sel_preset != "(nenhum)":
        try:
            preset = json.loads((presets_dir / sel_preset).read_text(encoding="utf-8"))
            for k, v in preset.items():
                st.session_state[k] = v
            st.success(f"Preset '{sel_preset}' carregado. Ajustes aplicados.")
        except Exception as e:
            st.error(f"Falha ao carregar preset: {e}")

    st.markdown("---")
    st.subheader("Saída")
    st.number_input("Largura de linha (visualização)", min_value=40, max_value=200, step=2, value=80, key="wrap_width")
    st.text_input("Pasta de artefatos (novo run)", value=base_runs_dir, key="artifacts_dir")

    st.markdown("---")
    st.subheader("Prompt Effatha — parâmetros")
    default_goal = st.session_state.get("last_prompt_goal", "Descrever impacto esperado e priorizar alvos para o tratamento especificado.")
    st.text_area("🎯 Objetivo do tratamento (texto livre)", value=default_goal, key="treatment_goal", height=120)
    mech = st.session_state.get("last_mech_flags", {"aproximar": True, "distanciar": True, "mimetizar": True})
    cma, cmd, cmm = st.columns(3)
    with cma:
        st.checkbox("Permitir Aproximar", value=bool(mech.get("aproximar", True)), key="allow_aproximar")
    with cmd:
        st.checkbox("Permitir Distanciar", value=bool(mech.get("distanciar", True)), key="allow_distanciar")
    with cmm:
        st.checkbox("Permitir Mimetizar", value=bool(mech.get("mimetizar", True)), key="allow_mimetizar")

    run_btn = st.button("🚀 Executar análise", type="primary", key="run_btn")

# -------------------------------
# Abas
# -------------------------------
(
    tab_resumo,
    tab_uniprot,
    tab_pdb,
    tab_ncbi_prot,
    tab_ncbi_gene,
    tab_regioes,
    tab_seqs,
    tab_blast_audit,
    tab_logs,
    tab_downloads
) = st.tabs(
    ["Resumo", "UniProt", "PDB", "NCBI Protein", "NCBI Gene", "Regiões", "Sequências", "Auditoria BLAST", "Logs", "Downloads"]
)

log_placeholder = tab_logs.empty()
summary_placeholder = tab_resumo.empty()
regions_placeholder = tab_regioes.empty()
uniprot_container = tab_uniprot.container()
pdb_container = tab_pdb.container()
ncbi_prot_container = tab_ncbi_prot.container()
ncbi_gene_container = tab_ncbi_gene.container()
seqs_container = tab_seqs.container()
blast_audit_container = tab_blast_audit.container()
downloads_container = tab_downloads.container()

# -------------------------------
# Render helpers
# -------------------------------
def kpi_source_box(container, context_like, title):
    with container:
        st.subheader(title)
        if not context_like:
            st.info("Sem dados.")
            return
        regs = context_like.get("regions", []) or []
        nt_segments = context_like.get("nt_segments", []) or []
        aa_any = context_like.get("aa_regions", []) or []
        dens_list = []
        for a in aa_any:
            d,_ = bracket_density(a.get("aa",""))
            dens_list.append(d)
        avg_d = (sum(dens_list)/len(dens_list)) if dens_list else 0.0
        st.markdown(
            f"""
            <div class="kpibox">
            <b>Regiões:</b> {len(regs)} &nbsp;|&nbsp; <b>Com NT:</b> {len(nt_segments)} &nbsp;|&nbsp; <b>AA com colchetes (média):</b> {avg_d:.1%}
            </div>
            """, unsafe_allow_html=True
        )

def render_density_plot(flags, title):
    if not plt or not flags:
        return
    vals = moving_avg(flags, k=15)
    fig, ax = plt.subplots()
    ax.plot(vals)
    ax.set_ylim(0, 1)
    ax.set_ylabel("densidade (média móvel)")
    ax.set_title(title)
    st.pyplot(fig, clear_figure=True)

# ===== BLOCO: Busca & Realce por Padrões =====
def patterns_ui_block():
    st.subheader("🔎 Busca & marcação por padrões (regex + catálogos)")
    cols = st.columns([1.5, 1.3, 1.2, 1.0])
    with cols[0]:
        sel_categories = st.multiselect("Categorias", options=list(PATTERN_CATALOG.keys()), default=["RNAi (siRNA/shRNA)"])
    with cols[1]:
        regex_manual = st.text_input("Regex manual (opcional)", value="")
    with cols[2]:
        mirna_seed = st.text_input("miRNA seed (5'→3', ex.: GAGGUAG)", value="")
    with cols[3]:
        mirna_scope = st.selectbox("miRNA alvo em", options=["mRNA (U)", "DNA (T)"], index=0)

    # Limites GC% (apenas anotação no tooltip; filtro visual é manual)
    cgc1, cgc2, _ = st.columns([1,1,1.2])
    with cgc1:
        gc_min = st.slider("GC% min (siRNA)", 0.0, 100.0, 30.0, 1.0)
    with cgc2:
        gc_max = st.slider("GC% max (siRNA)", 0.0, 100.0, 52.0, 1.0)

    # Monta lista de specs ativos
    active_specs = []
    for cat in sel_categories:
        active_specs.extend(PATTERN_CATALOG.get(cat, []))

    # Gera padrões para miRNA seed
    if mirna_seed.strip():
        seed = re.sub(r"[^ACGTUacgtu]", "", mirna_seed.strip())
        if seed:
            if mirna_scope.startswith("mRNA"):
                # alvo em mRNA — usamos complemento do seed (RNA)
                seed2_8 = seed[1:8] if len(seed) >= 8 else seed[1:]
                comp_2_8 = revcomp_rna(seed2_8)
                # 7mer-m8: match dos 7 nt (2–8)
                active_specs.append({
                    "key":"MIRNA_7m8_RNA",
                    "label":f"miRNA 7mer-m8 ({seed})",
                    "regex":comp_2_8,
                    "flags":re.I,
                    "scope":"mRNA",
                    "category_tag":"RNAi",
                    "tooltip":f"miRNA seed 7mer-m8 para {seed}",
                    "priority":12
                })
                # 7mer-A1: A upstream + 2–7
                if len(seed) >= 7:
                    comp_2_7 = revcomp_rna(seed[1:7])
                    active_specs.append({
                        "key":"MIRNA_7A1_RNA",
                        "label":f"miRNA 7mer-A1 ({seed})",
                        "regex":"A" + comp_2_7,
                        "flags":re.I,
                        "scope":"mRNA",
                        "category_tag":"RNAi",
                        "tooltip":f"miRNA seed 7mer-A1 (A1 + 2–7) para {seed}",
                        "priority":12
                    })
                # 8mer: A upstream + 2–8
                if len(seed) >= 8:
                    active_specs.append({
                        "key":"MIRNA_8MER_RNA",
                        "label":f"miRNA 8mer ({seed})",
                        "regex":"A" + comp_2_8,
                        "flags":re.I,
                        "scope":"mRNA",
                        "category_tag":"RNAi",
                        "tooltip":f"miRNA seed 8mer (A1 + 2–8) para {seed}",
                        "priority":12
                    })
            else:
                # alvo em DNA — complemento em DNA
                seed2_8 = seed[1:8] if len(seed) >= 8 else seed[1:]
                comp_2_8 = revcomp_dna(seed2_8)
                active_specs.append({
                    "key":"MIRNA_7m8_DNA",
                    "label":f"miRNA 7mer-m8 ({seed}) DNA",
                    "regex":comp_2_8,
                    "flags":re.I,
                    "scope":"DNA",
                    "category_tag":"RNAi",
                    "tooltip":f"miRNA seed 7mer-m8 (DNA) para {seed}",
                    "priority":12
                })
                if len(seed) >= 7:
                    comp_2_7 = revcomp_dna(seed[1:7])
                    active_specs.append({
                        "key":"MIRNA_7A1_DNA",
                        "label":f"miRNA 7mer-A1 ({seed}) DNA",
                        "regex":"A" + comp_2_7,
                        "flags":re.I,
                        "scope":"DNA",
                        "category_tag":"RNAi",
                        "tooltip":f"miRNA seed 7mer-A1 (DNA) para {seed}",
                        "priority":12
                    })
                if len(seed) >= 8:
                    active_specs.append({
                        "key":"MIRNA_8MER_DNA",
                        "label":f"miRNA 8mer ({seed}) DNA",
                        "regex":"A" + comp_2_8,
                        "flags":re.I,
                        "scope":"DNA",
                        "category_tag":"RNAi",
                        "tooltip":f"miRNA seed 8mer (DNA) para {seed}",
                        "priority":12
                    })

    # Regex manual (vira padrão "custom")
    if regex_manual.strip():
        active_specs.append({
            "key":"CUSTOM_REGEX",
            "label":"Regex manual",
            "regex":regex_manual.strip(),
            "flags":re.IGNORECASE,
            "scope":"AA",   # mostramos em todos; na hora do render filtramos por escopo
            "category_tag":"CUSTOM",
            "tooltip":"Busca manual (regex)",
            "priority":90
        })
        # Duplicamos para DNA/mRNA também
        active_specs.append({**active_specs[-1], "scope":"DNA"})
        active_specs.append({**active_specs[-1], "scope":"mRNA"})

    # Painel de seleção fina por padrão (checklist)
    if active_specs:
        st.caption("Seleção fina (habilite/desabilite padrões abaixo):")
        enabled = []
        for sp in active_specs:
            key = f"pat_{sp['key']}_{sp['scope']}"
            chk = st.checkbox(f"{sp['label']} — [{sp['scope']}]", value=True, key=key)
            if chk:
                enabled.append(sp)
        active_specs = enabled

    # Legenda de cores
    st.markdown(
        """
        **Cores:** <span class="hl-aa">AA</span> · <span class="hl-nt">NT</span> · <span class="hl-rnai">RNAi</span> ·
        <span class="hl-prom">Prom/Boxes</span> · <span class="hl-motif">Motivos AA</span> · <span class="hl-custom">Custom</span>
        """, unsafe_allow_html=True
    )
    st.caption("Passe o mouse sobre um trecho marcado para ver o tooltip (padrão/GC%/escopo).")

    return active_specs, gc_min, gc_max

def seq_block_markdown(label, seq_text, scope, specs_all, gc_min=0.0, gc_max=100.0):
    """Renderiza bloco com realce multicor e tooltips. 'scope' orienta o filtro de padrões."""
    if not seq_text:
        return
    st.markdown(f"**{label}**")
    # Filtra por escopo
    use_specs = []
    for sp in specs_all or []:
        sc = sp.get("scope")
        if scope == sc or (scope in ("DNA","mRNA") and sc=="NT"):
            # para siRNA: se tooltip guardar GC%, reforçar texto com alerta quando fora do range
            tooltip = sp.get("tooltip")
            if callable(tooltip) and "RNAi" in sp.get("category_tag",""):
                # embrulhamos o tooltip para acrescentar faixa GC quando necessário (apenas para os siRNA AA N19 TT/UU)
                def make_tt(tp):
                    def f(mtxt):
                        base = tp(mtxt)
                        # tenta extrair GC% do final
                        try:
                            # para AA[...]{19}TT pegamos mtxt[2:-2] ou correspondente
                            core = mtxt
                            if len(mtxt) >= 23 and (mtxt[:2].upper()=="AA") and (mtxt[-2:].upper() in ("TT","UU")):
                                core = mtxt[2:-2]
                            gcv = gc_percent(core)
                            flag = "" if (gc_min <= gcv <= gc_max) else "  ⚠ GC% fora do alvo"
                            return f"{base}{flag}"
                        except Exception:
                            return base
                    return f
                sp = {**sp, "tooltip": make_tt(tooltip)}
            use_specs.append(sp)
    compiled = compile_specs(use_specs)
    html = apply_multi_highlight(seq_text, compiled)
    st.markdown(html, unsafe_allow_html=True)

# ---- Mapeador proteína→códon
STD_CODON_TABLE = {
    'A': ['GCT','GCC','GCA','GCG'],
    'C': ['TGT','TGC'],
    'D': ['GAT','GAC'],
    'E': ['GAA','GAG'],
    'F': ['TTT','TTC'],
    'G': ['GGT','GGC','GGA','GGG'],
    'H': ['CAT','CAC'],
    'I': ['ATT','ATC','ATA'],
    'K': ['AAA','AAG'],
    'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
    'M': ['ATG'],
    'N': ['AAT','AAC'],
    'P': ['CCT','CCC','CCA','CCG'],
    'Q': ['CAA','CAG'],
    'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
    'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
    'T': ['ACT','ACC','ACA','ACG'],
    'V': ['GTT','GTC','GTA','GTG'],
    'W': ['TGG'],
    'Y': ['TAT','TAC'],
    '*': ['TAA','TAG','TGA']
}

def render_codon_mapper(context_like, container):
    if not context_like:
        with container:
            st.info("Sem NT disponível para mapear códons.")
        return
    nts = context_like.get("nt_segments", []) or []
    aa_any = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context_like.get("aa_regions", []) or []}
    if not nts:
        with container:
            st.info("Sem NT disponível para mapear códons.")
        return

    with container:
        st.subheader("Mapeador proteína → códon (por região)")
        options = [f"{s['tag']} {s['start']}-{s['end']}" for s in nts]
        sel = st.selectbox("Escolha a região", options=options, index=0, key="codon_region_sel")
        seg = nts[options.index(sel)]
        tag, start, end = seg["tag"], seg["start"], seg["end"]
        L = end - start + 1
        pos = st.number_input("Posição de AA dentro da região (1-based)", min_value=1, max_value=L, value=1, step=1, key="codon_pos_idx")
        dna = seg.get("dna","")
        mrna = seg.get("mrna","")
        aa = aa_any.get((tag, start, end), "")

        # pega AA ref na posição (considerando colchetes Effatha)
        idxAA = 0
        kpos = 1
        aa_ref = None
        while idxAA < len(aa) and kpos <= pos:
            ch = aa[idxAA]
            if ch == "[":
                j = aa.find("]", idxAA+1)
                token = aa[idxAA+1:j] if j != -1 else ""
                aa_ref = token.split("/")[0] if token else None
                idxAA = (j+1) if j != -1 else (idxAA+1)
            elif ch == "*" and idxAA == len(aa)-1:
                aa_ref = "*"; idxAA += 1
            else:
                aa_ref = ch
                idxAA += 1
            kpos += 1

        idx0 = (pos-1)*3
        codon_dna = dna[idx0:idx0+3] if len(dna) >= idx0+3 else ""
        codon_mrna = mrna[idx0:idx0+3] if len(mrna) >= idx0+3 else ""
        st.markdown(f"**AA(ref)**: `{aa_ref or '?'}` &nbsp; | &nbsp; **DNA(códon)**: `{codon_dna}` &nbsp; | &nbsp; **mRNA(códon)**: `{codon_mrna}`")

        aa_upper = (aa_ref or "").upper()
        syn = STD_CODON_TABLE.get(aa_upper, [])
        if syn:
            st.caption("Códons sinônimos (tabela padrão 1): " + ", ".join(syn))

# -------------------------------
# Render global
# -------------------------------
def render_region_expanders(context_like, container, title=None,
                            show_obs=True, show_fil=True, show_any=True,
                            only_favs=False,
                            active_specs=None, gc_min=0.0, gc_max=100.0):
    if context_like is None:
        with container:
            st.info("Dados não disponíveis para esta fonte.")
        return

    aa_idx_obs = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context_like.get("aa_regions_obs", []) or []}
    aa_idx_fil = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context_like.get("aa_regions_filtered", []) or []}
    aa_idx_any = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context_like.get("aa_regions", []) or []}
    cards = context_like.get("region_cards", []) or []
    nt_segments = context_like.get("nt_segments", []) or []

    rows = []
    if nt_segments:
        for seg in nt_segments:
            rows.append(("__NT__", seg["tag"], seg["start"], seg["end"], seg))
    else:
        for rec in context_like.get("aa_regions", []) or []:
            rows.append(("__AA__", rec["tag"], rec["start"], rec["end"], rec))

    if only_favs:
        rows = [r for r in rows if is_fav(r[1], r[2], r[3])]

    if title:
        container.subheader(title)

    if not rows:
        container.info("Sem sequências para mostrar.")
        return

    # Controles locais da seção
    c1, c2, c3 = container.columns([1.2,1,1])
    with c1:
        exp_all = st.checkbox("Expandir todos", value=st.session_state.get("expand_all", False), key=f"expand_all_{title or 'main'}")
        st.session_state["expand_all"] = exp_all
    with c2:
        favs_only_local = st.checkbox("⭐ Somente favoritas", value=only_favs, key=f"only_favs_{title or 'main'}")
    with c3:
        st.caption("Use a seção acima para escolher padrões/regex.")

    for kind, tag, start, end, payload in rows:
        c = pick_card(cards, tag, start, end)
        header = f"**{tag}** {start}-{end}"
        if c and c.get("length"):
            header += f" (len={c.get('length')})"
        if c and c.get("note"):
            header += f" — {c.get('note')}"
        badges = []
        if tag.upper() == "FULL":
            badges.append('<span class="badge full">FULL</span>')
        if c and ("uniprot" in (c.get("note","").lower())):
            badges.append('<span class="badge uniprot">UniProt</span>')
        if c and ("genpept" in (c.get("note","").lower())):
            badges.append('<span class="badge genpept">GenPept</span>')
        if c and c.get("ligand_name"):
            badges.append('<span class="badge lig">Lig.</span>')

        with container.expander(f"{tag} {start}-{end}", expanded=st.session_state.get("expand_all", False)):
            st.markdown("".join(badges) + " " + header, unsafe_allow_html=True)
            fav_state = is_fav(tag, start, end)
            fav_toggle = st.checkbox("⭐ Favoritar", value=fav_state, key=f"fav_{title}_{tag}_{start}_{end}")
            toggle_fav(tag, start, end, fav_toggle)

            # AA
            aa_obs = aa_idx_obs.get((tag, start, end), "")
            aa_fil = aa_idx_fil.get((tag, start, end), "")
            aa_any = "" if (aa_obs or aa_fil) else aa_idx_any.get((tag, start, end), "")

            if show_obs and aa_obs:
                seq_block_markdown("AA (observado)", aa_obs, scope="AA", specs_all=active_specs, gc_min=gc_min, gc_max=gc_max)
                d, flags = bracket_density(aa_obs)
                st.caption(f"Densidade de colchetes (obs): {d:.1%}")
                render_density_plot(flags, f"{tag} {start}-{end} — AA(obs)")

            if show_fil and aa_fil:
                seq_block_markdown("AA (filtrado)", aa_fil, scope="AA", specs_all=active_specs, gc_min=gc_min, gc_max=gc_max)
                d, flags = bracket_density(aa_fil)
                st.caption(f"Densidade de colchetes (filtr.): {d:.1%}")
                render_density_plot(flags, f"{tag} {start}-{end} — AA(filt)")

            if show_any and aa_any:
                seq_block_markdown("AA", aa_any, scope="AA", specs_all=active_specs, gc_min=gc_min, gc_max=gc_max)
                d, flags = bracket_density(aa_any)
                st.caption(f"Densidade de colchetes: {d:.1%}")
                render_density_plot(flags, f"{tag} {start}-{end} — AA")

            if kind == "__NT__":
                dna = payload.get("dna","")
                mrna = payload.get("mrna","")
                seq_block_markdown("DNA", dna, scope="DNA", specs_all=active_specs, gc_min=gc_min, gc_max=gc_max)
                dD, fD = bracket_density(dna)
                st.caption(f"Densidade de colchetes (DNA): {dD:.1%}")
                render_density_plot(fD, f"{tag} {start}-{end} — DNA")

                seq_block_markdown("mRNA", mrna, scope="mRNA", specs_all=active_specs, gc_min=gc_min, gc_max=gc_max)
                dR, fR = bracket_density(mrna)
                st.caption(f"Densidade de colchetes (mRNA): {dR:.1%}")
                render_density_plot(fR, f"{tag} {start}-{end} — mRNA")

def render_outputs(context, paths, enable_blast_flag):
    protein = context.get("protein", {}) or {}
    layers = context.get("layers", {}) or {}
    blast = context.get("blast", {}) or {}

    with summary_placeholder:
        st.subheader("Resumo da Execução")
        colA, colB, colC = st.columns(3)
        with colA:
            st.metric("Acesso", protein.get("accession","?"))
            st.metric("Tamanho (aa)", protein.get("length","?"))
        with colB:
            st.metric("Espécie", f"{protein.get('organism','?')} (taxid={protein.get('taxid','?')})")
            links = build_quick_links(protein)
            if links:
                st.markdown("Links: " + " | ".join([f"[{n}]({u})" for n,u in links]))
        with colC:
            st.metric("Regiões-alvo", len(context.get("regions", [])))

        st.markdown("**Camadas ativas**")
        st.caption(
            f"UniProt VARIANT = {'ON' if layers.get('uniprot_variant_enabled') else 'OFF'} | "
            f"Proteins Variation = {'ON' if layers.get('proteins_variation_enabled') else 'OFF'} | "
            f"BLAST = {'ON' if layers.get('blast_enabled') else 'OFF'}"
        )

        st.markdown("**Variações — contagens**")
        c1, c2, c3, c4 = st.columns(4)
        with c1: st.metric("BLAST (pos.)", blast.get("variant_positions_blast", 0))
        with c2: st.metric("UniProt VARIANT (pos.)", blast.get("variant_positions_uniprot", 0))
        with c3: st.metric("Proteins Variation (pos.)", blast.get("variant_positions_proteins", 0))
        with c4: st.metric("FILTRADAS (pos.)", blast.get("variant_positions_filtered", 0))

    srcs = context.get("sources", {}) or {}
    kpi_source_box(uniprot_container, srcs.get("uniprot"), "UniProt")
    kpi_source_box(pdb_container, srcs.get("pdb"), "PDB")
    kpi_source_box(ncbi_prot_container, srcs.get("ncbi_protein"), "NCBI Protein")
    kpi_source_box(ncbi_gene_container, srcs.get("ncbi_gene"), "NCBI Gene")

    # ----- Regiões (cards)
    cards = context.get("region_cards", []) or []
    with regions_placeholder:
        if not cards or pd is None:
            st.info("Sem region cards disponíveis.")
        else:
            rows = []
            for c in cards:
                rows.append({
                    "Tag": c.get("tag"),
                    "Início": c.get("start"),
                    "Fim": c.get("end"),
                    "len": c.get("length"),
                    "Descrição": (c.get("note") or ""),
                    "Ligante": (c.get("ligand_name") or ""),
                    "ChEBI": (c.get("ligand_chebi") or ""),
                    "Íon metálico": ("sim" if c.get("ligand_is_metal") else ("não" if c.get("ligand_is_metal") is not None else "")),
                    "% var (obs)": f"{100.0*(c.get('pct_var_obs') or 0):.1f}%",
                    "% var (filtrado)": f"{100.0*(c.get('pct_var_fil') or 0):.1f}%",
                    "AF máx": (f"{c.get('max_af'):.4f}" if c.get("max_af") is not None else ""),
                    "Clínica": ", ".join([f"{k}={v}" for k,v in (c.get("clin_counts") or {}).items()]),
                    "Evidência manual (pos.)": c.get("pos_with_manual_ev"),
                    "SIFT (tol/del)": f"{c.get('sift_tol')}/{c.get('sift_del')}",
                    "PolyPhen (ben/dam)": f"{c.get('poly_benign')}/{c.get('poly_damaging')}",
                    "⭐": "★" if is_fav(c.get("tag"), c.get("start"), c.get("end")) else ""
                })
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True)

    # ===== Busca & marcação (novo bloco)
    with seqs_container:
        active_specs, gc_min, gc_max = patterns_ui_block()

        st.markdown("---")
        st.subheader("Sequências (vista consolidada)")
        cobs, cfil, candy = st.columns(3)
        with cobs:
            show_obs = st.checkbox("Mostrar AA (observado)", value=True, key="show_obs_main")
        with cfil:
            show_fil = st.checkbox("Mostrar AA (filtrado)", value=True, key="show_fil_main")
        with candy:
            show_any = st.checkbox("Mostrar AA (bruto)", value=True, key="show_any_main")

        render_region_expanders(
            context, st.container(), title="Consolidado",
            show_obs=show_obs, show_fil=show_fil, show_any=show_any,
            only_favs=False,
            active_specs=active_specs, gc_min=gc_min, gc_max=gc_max
        )

        st.markdown("---")
        st.subheader("Sequências por fonte")
        col_f1, col_f2 = st.columns(2)
        with col_f1:
            st.caption("UniProt")
            render_region_expanders(srcs.get("uniprot"), st.container(), title="UniProt",
                                    show_obs=show_obs, show_fil=show_fil, show_any=show_any,
                                    only_favs=False, active_specs=active_specs, gc_min=gc_min, gc_max=gc_max)
            st.caption("NCBI Protein")
            render_region_expanders(srcs.get("ncbi_protein"), st.container(), title="NCBI Protein",
                                    show_obs=show_obs, show_fil=show_fil, show_any=show_any,
                                    only_favs=False, active_specs=active_specs, gc_min=gc_min, gc_max=gc_max)
        with col_f2:
            st.caption("PDB")
            render_region_expanders(srcs.get("pdb"), st.container(), title="PDB",
                                    show_obs=show_obs, show_fil=show_fil, show_any=show_any,
                                    only_favs=False, active_specs=active_specs, gc_min=gc_min, gc_max=gc_max)
            st.caption("NCBI Gene")
            render_region_expanders(srcs.get("ncbi_gene"), st.container(), title="NCBI Gene",
                                    show_obs=show_obs, show_fil=show_fil, show_any=show_any,
                                    only_favs=False, active_specs=active_specs, gc_min=gc_min, gc_max=gc_max)

        st.markdown("---")
        st.subheader("Ferramentas")
        tool_col1, tool_col2 = st.columns(2)
        with tool_col1:
            render_codon_mapper(context, st.container())
        with tool_col2:
            st.info("Dica: combine categorias (ex.: RNAi + Promotores) e regex manual para investigações específicas. Passe o mouse nos trechos marcados para ver o tooltip.")

    # ----- Auditoria BLAST
    with blast_audit_container:
        st.subheader("Auditoria BLAST — hits aceitos")
        run_dir = Path(st.session_state.get("last_run_dir") or "")
        if not run_dir.exists():
            st.info("Execute um run ou carregue um do histórico para ver auditoria.")
        else:
            csvs = list(run_dir.glob("blast_hits*.csv")) + list(run_dir.glob("blast_hits_*.csv"))
            if not csvs:
                st.info("Nenhum arquivo blast_hits*.csv encontrado neste run.")
            else:
                if pd is None:
                    st.warning("Pandas não disponível para exibir auditoria.")
                else:
                    dfs = []
                    for p in csvs:
                        try:
                            df = pd.read_csv(p)
                            df["arquivo"] = p.name
                            dfs.append(df)
                        except Exception:
                            pass
                    if not dfs:
                        st.info("Não foi possível ler os CSVs de auditoria BLAST.")
                    else:
                        df_all = pd.concat(dfs, ignore_index=True)
                        # Normaliza nomes comuns
                        rename_map = {
                            "region":"region_tag",
                            "percent_identity":"pct_identity",
                            "sbjct_id":"accession"  # se vier outro nome
                        }
                        for c_old, c_new in rename_map.items():
                            if c_old in df_all.columns and c_new not in df_all.columns:
                                df_all[c_new] = df_all[c_old]

                        cols = df_all.columns.tolist()
                        def safe(name): return name if name in cols else None
                        col_region = safe("region_tag")
                        col_kind = safe("kind")
                        col_db = safe("db")
                        col_pid = safe("pct_identity")
                        col_qcov = safe("query_coverage")
                        col_acc = safe("accession")

                        cA, cB, cC, cD = st.columns(4)
                        with cA:
                            region_sel = st.multiselect("Região", sorted(df_all[col_region].dropna().unique().tolist()) if col_region else [], default=None)
                        with cB:
                            kind_sel = st.multiselect("Nível (AA/NT)", sorted(df_all[col_kind].dropna().unique().tolist()) if col_kind else [], default=None)
                        with cC:
                            min_id = st.slider("Min %ID", 0.0, 100.0, 90.0, 0.5)
                        with cD:
                            min_qc = st.slider("Min %COV", 0.0, 100.0, 90.0, 0.5)

                        dfv = df_all.copy()
                        if col_region and region_sel:
                            dfv = dfv[dfv[col_region].isin(region_sel)]
                        if col_kind and kind_sel:
                            dfv = dfv[dfv[col_kind].isin(kind_sel)]
                        if col_pid:
                            dfv = dfv[dfv[col_pid] >= min_id]
                        if col_qcov:
                            dfv = dfv[dfv[col_qcov] >= min_qc]

                        if col_acc:
                            links = []
                            for s in dfv[col_acc].astype(str).tolist():
                                kind, acc = guess_accession_from_sbjct(s)
                                links.append(make_ext_link(kind, acc))
                            dfv = dfv.assign(link=links)

                        st.dataframe(dfv, use_container_width=True)
                        st.download_button(
                            "Exportar auditoria filtrada (.csv)",
                            data=dfv.to_csv(index=False).encode("utf-8"),
                            file_name="blast_audit_filtrada.csv",
                            mime="text/csv"
                        )

    # ----- Downloads
    with downloads_container:
        st.subheader("Downloads")
        try:
            rp = Path(paths.get("report","")) if paths else None
            if rp and rp.exists():
                st.download_button("Baixar relatório (.txt)", data=rp.read_bytes(),
                                   file_name=rp.name, mime="text/plain")

            rcsv = Path(paths.get("regions_csv","")) if paths else None
            if rcsv and rcsv.exists():
                st.download_button("Baixar regiões (.csv)", data=rcsv.read_bytes(),
                                   file_name=rcsv.name, mime="text/csv")

            bcsv = Path(paths.get("blast_csv","")) if paths else None
            if bcsv and bcsv.exists() and enable_blast_flag:
                st.download_button("Baixar variantes BLAST (.csv)", data=bcsv.read_bytes(),
                                   file_name=bcsv.name, mime="text/csv")

            ctx_path = Path(st.session_state.get("last_run_dir","")) / "context_summary.json" if st.session_state.get("last_run_dir") else None
            if ctx_path and ctx_path.exists():
                st.download_button("Baixar contexto (.json)", data=ctx_path.read_bytes(),
                                   file_name=ctx_path.name, mime="application/json")

            run_dir = Path(st.session_state.get("last_run_dir") or "")
            if run_dir.exists():
                zip_path = run_dir.with_suffix(".zip")
                with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
                    for path in run_dir.rglob("*"):
                        z.write(path, arcname=path.relative_to(run_dir))
                st.download_button("Baixar bundle do run (.zip)", data=zip_path.read_bytes(),
                                   file_name=zip_path.name, mime="application/zip")
        except Exception as e:
            st.error(f"Erro ao preparar downloads: {e}")

# -------------------------------
# Prompt Effatha (arquivo) — igual ao anterior (resumo)
# -------------------------------
def _list_regions_lines(context):
    lines = []
    for r in context.get("regions", []) or []:
        note = f" — {r.get('note','')}" if r.get('note') else ""
        src = r.get('source'); typ = r.get('type'); tag = r.get('tag')
        if src or typ:
            lines.append(f"- {src or ''}:{typ or ''} {r.get('start')}-{r.get('end')}{note}")
        else:
            lines.append(f"- {tag or ''} {r.get('start')}-{r.get('end')}{note}")
    return "\n".join(lines) if lines else "(sem regiões)"

def _seq_block(label, seq_text):
    if not seq_text:
        return ""
    return f"**{label}**\n\n```\n{seq_text}\n```\n"

def _render_region_sequences(context):
    aa_obs = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions_obs", [])}
    aa_fil = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions_filtered", [])}
    aa_any = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions", [])}

    blocks = []
    for seg in context.get("nt_segments", []) or []:
        tag = seg["tag"]; start = seg["start"]; end = seg["end"]
        note = seg.get("note") or ""
        header = f"### {tag} {start}-{end}"
        if note:
            header += f" — {note}"
        aa1 = aa_obs.get((tag,start,end), "")
        aa2 = aa_fil.get((tag,start,end), "")
        aa3 = "" if (aa1 or aa2) else aa_any.get((tag,start,end), "")
        dna = seg.get("dna","")
        mrna = seg.get("mrna","")
        section = [header, ""]
        section.append(_seq_block("AA (observado)", aa1))
        section.append(_seq_block("AA (filtrado)", aa2))
        section.append(_seq_block("AA", aa3))
        section.append(_seq_block("DNA", dna))
        section.append(_seq_block("mRNA", mrna))
        blocks.append("\n".join([x for x in section if x]))
    return "\n\n".join(blocks) if blocks else "_Sem sequências por região disponíveis._"

def write_prompt_file(run_dir, *, context, treatment_goal, mech_flags, report_path, slug, ts):
    protein = context.get("protein", {}) or {}
    blast = context.get("blast", {}) or {}
    layers = context.get("layers", {}) or {}

    md = f"""# Prompt Effatha — Análise e Priorização de Sequências (para GPT Pro)

**Objetivo do tratamento:** {treatment_goal or "(não informado)"}  
**Mecanismos permitidos nesta análise:** {", ".join([k for k,v in mech_flags.items() if v]) or "aproximar, distanciar e mimetizar"}

## 0) Papel da IA
Você é uma IA especialista em biologia molecular/sistemas, treinada para avaliar **potencial de intervenção por frequências Effatha** em proteínas/genes.

Para **cada sequência** abaixo:
1. Classifique o **potencial** Effatha como **"Alto"**, **"Moderado"**, **"Baixo"** ou **"Sem potencial"**.
2. Indique **qual mecanismo** é mais promissor (**aproximar**, **distanciar**, **mimetizar**) — restrito aos **permitidos** acima.
3. Explique em **≤ 8 linhas**, citando região/função, contexto estrutural, **variantes por posição** (`[ref/alt/…]`), efeitos esperados e riscos.
4. Entregue um **ranking global** (Top N) e uma **lista de descarte** (“Sem potencial”) com 1 linha de motivo.
5. Use **apenas dados reais**; se faltar, retorne **"dado indisponível"** (não inventar).

### Regras de Sintaxe e Processo
- **Sintaxe Effatha** por **posição** (AA **e** NT): `[A/L/E]` — **NT é por base, não por trinca**.
- O Analyzer:
  - Expande regiões com **flancos** (ex.: ±5 aa).
  - Insere **variantes observadas** (BLAST/UniProt/Proteins Variation) **por posição** no formato Effatha.
  - Inclui **DNA/mRNA** somente quando há **CDS real mapeado** (sem retrotradução/heurística).

## 1) Contexto do caso
- **Acesso**: {protein.get('accession','?')}
- **Organismo**: {protein.get('organism','?')} (taxid={protein.get('taxid','?')})
- **Tamanho (aa)**: {protein.get('length','?')}
- **Camadas ativas no Analyzer**:
  - UniProt VARIANT = { "ON" if layers.get("uniprot_variant_enabled") else "OFF" }
  - Proteins Variation API = { "ON" if layers.get("proteins_variation_enabled") else "OFF" }
  - BLAST = { "ON" if layers.get("blast_enabled") else "OFF" }
- **BLAST (resumo)**: pos_var(blast)={blast.get('variant_positions_blast',0)}

## 2) Regiões-alvo (após expansão)
{_list_regions_lines(context)}

## 3) Sequências por região (AA + DNA/mRNA)
{_render_region_sequences(context)}
"""
    prompt_md_path = Path(run_dir) / f"prompt_effatha_{slug}_{ts}.md"
    prompt_txt_path = Path(run_dir) / f"prompt_effatha_{slug}_{ts}.txt"
    prompt_md_path.write_text(md, encoding="utf-8")
    prompt_txt_path.write_text(md, encoding="utf-8")
    return str(prompt_md_path), str(prompt_txt_path)

# -------------------------------
# Execução
# -------------------------------
def _sanitize(s: str) -> str:
    return "".join(ch for ch in (s or "") if ch.isalnum() or ch in ("-", "_", ".")).strip() or "run"

def save_session(run_dir, context, paths, enable_blast, prompt_goal, mech_flags):
    st.session_state["last_run_dir"] = str(run_dir)
    st.session_state["last_context"] = context
    st.session_state["last_paths"] = paths
    st.session_state["last_enable_blast"] = bool(enable_blast)
    st.session_state["last_prompt_goal"] = prompt_goal
    st.session_state["last_mech_flags"] = mech_flags

def run_pipeline_and_render():
    pp = Path(st.session_state["pipeline_path"])
    if not pp.is_file():
        st.error(f"Arquivo não encontrado: {pp}")
        st.stop()

    run_dir = make_run_dir(st.session_state.get("artifacts_dir") or "runs")
    src_label = st.session_state["source_choice"]
    if src_label == "UniProt":
        id_hint = st.session_state.get("uniprot_acc", "")
    elif src_label == "PDB":
        id_hint = st.session_state.get("pdb_id", "")
    elif src_label == "NCBI Protein":
        id_hint = st.session_state.get("ncbi_protein_acc", "")
    else:
        if st.session_state.get("gene_mode") == "GeneID":
            id_hint = st.session_state.get("gene_id", "")
        else:
            sym = st.session_state.get("gene_symbol", "")
            tax = str(st.session_state.get("gene_taxid") or "")
            id_hint = f"{sym}_{tax}" if sym else tax

    slug = _sanitize(id_hint)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    report_path = str(run_dir / f"report_{slug}_{ts}.txt")
    regions_csv_path = str(run_dir / f"regions_{slug}_{ts}.csv")
    blast_csv_path = str(run_dir / f"variants_blast_{slug}_{ts}.csv")

    cfg_overrides = {
        "input": {
            "source": ("uniprot" if src_label=="UniProt" else
                       "pdb" if src_label=="PDB" else
                       "ncbi_protein" if src_label=="NCBI Protein" else
                       "ncbi_gene"),
            "uniprot_acc": st.session_state.get("uniprot_acc","").strip(),
            "pdb_id": st.session_state.get("pdb_id","").strip(),
            "ncbi_protein_acc": st.session_state.get("ncbi_protein_acc","").strip(),
            "gene": {
                "id_type": ("entrez" if st.session_state.get("gene_mode")=="GeneID" else "symbol"),
                "id": st.session_state.get("gene_id","").strip(),
                "symbol": st.session_state.get("gene_symbol","").strip(),
                "taxid": int(st.session_state.get("gene_taxid") or 9606),
                "isoform_policy": st.session_state.get("isoform_policy","longest"),
            },
        },
        "external_layers": {
            "uniprot_variant": {"enable": bool(st.session_state.get("enable_uniprot_variant", True))},
            "variation_api": {"enable": bool(st.session_state.get("enable_variation_api", True))},
        },
        "blast": {
            "enable": bool(st.session_state.get("enable_blast", False)),
            "protein": {
                "db": st.session_state.get("blast_db", "swissprot"),
                "hitlist_size": int(st.session_state.get("blast_hitlist_size", 25)),
                "expect": 1e-5,
                "min_identity": float(st.session_state.get("blast_min_id", 0.97)),
                "min_query_coverage": float(st.session_state.get("blast_min_cov", 0.95)),
            },
        },
        "variation_filters": {
            "enabled": bool(st.session_state.get("enable_filters", True)),
            "min_af": float(st.session_state.get("min_af", 0.01)),
            "clinical": {
                "allow": st.session_state.get("allow_clin", []),
                "allow_if_conflicting": bool(st.session_state.get("allow_if_conflicting", False)),
                "allow_uncertain": bool(st.session_state.get("allow_uncertain", False)),
            },
            "evidence": {
                "assertion": st.session_state.get("evidence_mode", "prefer_manual"),
                "treat_missing_as": st.session_state.get("missing_evid", "include"),
            },
            "predictions": {
                "sift_allow": [str(x) for x in st.session_state.get("sift_allow", [])],
                "polyphen_allow": [str(x) for x in st.session_state.get("polyphen_allow", [])],
                "require_both": bool(st.session_state.get("require_both", False)),
            },
            "include_blast_in_filtered": bool(st.session_state.get("include_blast_in_filtered", False)),
        },
        "regions": {
            "point_flank": int(st.session_state.get("point_flank", 5)),
            "default_min_len": int(st.session_state.get("default_min_len", 6)),
            "merge_overlaps": bool(st.session_state.get("merge_overlaps", True)),
            "add_full": bool(st.session_state.get("add_full", True)),
            "include_ncbi_protein_features": bool(st.session_state.get("include_ncbi_features", True)),
        },
        "output": {
            "report_txt": report_path,
            "export_regions_csv": regions_csv_path,
            "export_csv": (blast_csv_path if st.session_state.get("enable_blast", False) else None),
            "artifacts_dir": str(run_dir),
            "blast_progress": True,
        },
    }

    # opcional: restringir BLAST pela espécie (se o pipeline suportar)
    cfg_overrides["blast"]["same_species_only"] = bool(st.session_state.get("blast_same_species_only", True))

    # escrever shim e executar
    shim_path = write_shim(run_dir, st.session_state["pipeline_path"], cfg_overrides)
    py_exe = sys.executable

    with st.status("Executando pipeline…", expanded=True) as status:
        rc, out, err = run_subprocess(py_exe, shim_path, run_dir)
        log_placeholder.code((out or "") + ("\n" + err if err else ""))
        if rc == 0:
            status.update(label="Execução concluída", state="complete")
        else:
            status.update(label=f"Erro (exit code {rc}) — veja 'Logs'", state="error")

    # carregar contexto
    context_path = Path(run_dir) / "context_summary.json"
    context = None
    if context_path.exists():
        try:
            context = json.loads(context_path.read_text(encoding="utf-8"))
        except Exception as e:
            tab_resumo.error(f"Falha ao ler context_summary.json: {e}")

    # gerar prompt Effatha (arquivo)
    treatment_goal = st.session_state.get("treatment_goal", "").strip()
    mech_flags = {
        "aproximar": bool(st.session_state.get("allow_aproximar", True)),
        "distanciar": bool(st.session_state.get("allow_distanciar", True)),
        "mimetizar": bool(st.session_state.get("allow_mimetizar", True)),
    }
    prompt_md_path = ""
    prompt_txt_path = ""
    if context:
        try:
            prompt_md_path, prompt_txt_path = write_prompt_file(
                run_dir,
                context=context,
                treatment_goal=treatment_goal,
                mech_flags=mech_flags,
                report_path=report_path,
                slug=slug,
                ts=ts,
            )
        except Exception as e:
            tab_resumo.warning(f"Falha ao gerar Prompt Effatha: {e}")

    paths = {
        "report": report_path,
        "regions_csv": regions_csv_path,
        "blast_csv": blast_csv_path,
        "prompt_md": (prompt_md_path or ""),
        "prompt_txt": (prompt_txt_path or ""),
    }
    save_session(run_dir, context, paths, st.session_state.get("enable_blast", False), treatment_goal, mech_flags)

    if context:
        render_outputs(context, paths, st.session_state.get("last_enable_blast", False))
    else:
        tab_resumo.info("Contexto não encontrado — verifique os Logs para erros no pipeline.")

    # -------------------------------
    # Run / Restore
    # -------------------------------


if st.session_state.get("run_btn"):
    run_pipeline_and_render()
elif st.session_state.get("last_context"):
    render_outputs(
        st.session_state.get("last_context"),
        st.session_state.get("last_paths", {}),
        st.session_state.get("last_enable_blast", False),
    )
else:
    summary_placeholder.info("Pronto para executar. Configure os parâmetros e clique em **Executar análise**.")

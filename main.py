#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Effatha — Gene & Protein Analyzer (com regiões funcionais + gene gating)
- Normaliza entrada (UniProt, PDB, NCBI Protein, NCBI Gene) para um alvo UniProt
- Extrai features reais do UniProt (com flancos em sites pontuais)
- Mapeia CDS real (GenBank/nuccore) sem retrotradução
- BLASTp/BLASTn por FULL + por região → sintaxe Effatha [A/L/E] (AA) e [A/C/G/T] (NT)
- Filtros: mesma espécie (txid…[ORGN]) e **mesmo gene** (regex sobre título do hit, com *relax* por janela)
- Saídas: runs/report.txt, runs/regions.csv, runs/context_summary.json (+ opcional variants_blast.csv; auditoria runs/blast_hits.csv)

Requisitos:
  pip install biopython requests psutil
  export NCBI_EMAIL="seu_email@dominio"
  export NCBI_API_KEY="sua_chave"   # opcional, mas recomendado

Encerramento limpo:
- SIGINT/SIGTERM, atexit, fechamento de sessão HTTP e kill de filhos (psutil)
- HARD_EXIT=1 força os._exit
- EFFATHA_TIMEOUT_SECONDS=<seg> watchdog
"""

from __future__ import annotations

import os, io, re, csv, json, time, atexit, signal, sys, threading
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any, Set

import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

# ============== Encerramento / sinais / watchdog ==============

try:
    import psutil  # type: ignore
    _HAVE_PSUTIL = True
except Exception:
    _HAVE_PSUTIL = False

_HARD_EXIT = os.getenv("HARD_EXIT", "0") == "1"
_TIMEOUT_SECONDS = int(os.getenv("EFFATHA_TIMEOUT_SECONDS", "0") or "0")
_watchdog_timer: Optional[threading.Timer] = None

def _kill_children():
    if not _HAVE_PSUTIL:
        return
    try:
        me = psutil.Process(os.getpid())
        children = me.children(recursive=True)
        if children:
            print(f"[shim] Encerrando {len(children)} processo(s)-filho…", flush=True)
        for c in children:
            try:
                c.terminate()
            except Exception:
                pass
        gone, alive = psutil.wait_procs(children, timeout=5)
        for a in alive:
            try:
                a.kill()
            except Exception:
                pass
    except Exception:
        pass

def _final_flush():
    for s in (sys.stdout, sys.stderr):
        try: s.flush()
        except Exception: pass

def _finalize(code: int = 0):
    try: SESSION.close()
    except Exception: pass
    _kill_children()
    _final_flush()

def _finalize_and_exit(code: int = 0):
    _finalize(code)
    if _HARD_EXIT:
        os._exit(code)
    sys.exit(code)

def _signal_handler(signum, frame):
    print(f"[shim] Sinal {signum} — finalizando…", flush=True)
    _finalize_and_exit(130)

def _install_signal_handlers():
    for sig in (signal.SIGINT, signal.SIGTERM):
        try: signal.signal(sig, _signal_handler)
        except Exception: pass

def _start_watchdog():
    global _watchdog_timer
    if _TIMEOUT_SECONDS > 0:
        def _timeout():
            print(f"[shim] Tempo limite atingido ({_TIMEOUT_SECONDS}s). Encerrando.", flush=True)
            _finalize_and_exit(124)
        _watchdog_timer = threading.Timer(_TIMEOUT_SECONDS, _timeout)
        _watchdog_timer.daemon = True
        _watchdog_timer.start()

def _cancel_watchdog():
    global _watchdog_timer
    try:
        if _watchdog_timer: _watchdog_timer.cancel()
    except Exception:
        pass

@atexit.register
def _on_atexit():
    try: SESSION.close()
    except Exception: pass
    _kill_children()
    _final_flush()

# ========================= Config =========================
CONFIG: Dict[str, Any] = {
    "input": {
        "source": "ncbi_protein",  # "uniprot" | "pdb" | "ncbi_protein" | "ncbi_gene"
        "uniprot_acc": "P00533",
        "pdb_id": "5XNL",
        "ncbi_protein_acc": "AML61188.1",
        "gene": {"id_type": "entrez", "id": "1956", "symbol": "", "taxid": 9606, "isoform_policy": "longest"},
    },

    "regions": {
        "use_uniprot_features": True,
        "include_feature_types": [
            "domain","region","region of interest","repeat","coiled coil","zinc finger","motif","compositional bias",
            "transmembrane","topological domain","intramembrane",
            "active site","site","binding site","metal binding","calcium binding","dna binding","nucleotide binding",
            "signal peptide","transit peptide","propeptide","peptide","chain","initiator methionine",
            "glycosylation site","lipidation","modified residue","disulfide bond","cross-link",
            "natural variant","mutagenesis","sequence conflict","non-standard residue","non-terminal residue","non-adjacent residues",
            "helix","beta strand","turn"
        ],
        "point_flank": 5,
        "default_min_len": 6,
        "merge_overlaps": True,
        "add_full": True,
        "fallback_genpept_if_uniprot_featureless": True
    },

    "blast": {
        "enable": True,
        "same_species_only": True,
        "log_offspecies_preview": False,

        # ===== Protein (AA) =====
        "protein": {
            "dbs": ["refseq_protein", "nr"],
            "db": "refseq_protein",
            "hitlist_size": 200,
            "expect": 1e-5,
            "min_identity": 0.90,
            "min_query_coverage": 0.90,
            "entrez_query": None,

            # FULL→primeiro (semeia variantes)
            "smart_full_first": True,
            "full_dbs": None,
            "full_hitlist_size": 50,
            "full_min_query_coverage": 0.60,
            "max_alignments_process": 200,

            "regions_min_len_blast": 40,
            "gap_frac_threshold": 0.15,
            "skip_tags_for_region_blast": ["COMPOSITIONAL BIAS"],

            # baixa complexidade
            "filter_low_complexity": True,
        },

        # ===== Nucleotide (DNA/mRNA) =====
        "nt": {
            "dna_db": "nt",
            "rna_db": "refseq_rna",
            "hitlist_size": 25,
            "expect": 1e-10,
            "min_identity": 0.98,
            "min_query_coverage": 0.95,
            "megablast": True,
            "entrez_query": None,
            "filter_low_complexity": True,
        },

        # ===== Gene gating =====
        "gene_gate": {
            "enable": True,
            "accept_from_uniprot": True,   # usa símbolo/sinônimos do UniProt
            "extra_accept": [],            # strings adicionais (ex.: ["TP53"])
            "exclude": [],                 # strings para excluir (parálogos), ex.: ["CREBBP","CBP"]
            "relax_if_no_hits": True       # se nenhum hit passar, reprocessa janela ignorando gating
        }
    },

    "cds_mapping": {"enable": True},

    "output": {
        "artifacts_dir": "runs",
        "report_txt": "runs/report.txt",
        "export_regions_csv": "runs/regions.csv",
        "export_csv": None,  # "runs/variants_blast.csv"
        "blast_progress": True,
        "log_hits": True,
        "blast_hits_csv": "runs/blast_hits.csv",
    },

    "ncbi": {"email": "lovato.rodrigo@gmail.com", "api_key": ""}
}

# ========================= NCBI / HTTP =========================
NCBI_EMAIL   = os.getenv("NCBI_EMAIL", "").strip()
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()
NCBI_EMAIL   = CONFIG.get("ncbi", {}).get("email", NCBI_EMAIL).strip()
NCBI_API_KEY = CONFIG.get("ncbi", {}).get("api_key", NCBI_API_KEY).strip()

if NCBI_EMAIL: Entrez.email = NCBI_EMAIL
if NCBI_API_KEY: Entrez.api_key = NCBI_API_KEY
if not (NCBI_EMAIL or getattr(Entrez, "email", None)): Entrez.tool = "Effatha-GPA"
if not getattr(Entrez, "email", None):
    raise RuntimeError("NCBI_EMAIL não definido. Defina ENV NCBI_EMAIL ou CONFIG['ncbi']['email'].")

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Effatha-GPA/3.1 (+https://effatha)"})

def _http_retry(method: str, url: str, **kwargs):
    max_tries = kwargs.pop("max_tries", 6)
    backoff = 1.5
    for i in range(max_tries):
        try:
            resp = SESSION.request(method, url, timeout=45, **kwargs)
            if resp.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{resp.status_code} transient", response=resp)
            resp.raise_for_status()
            return resp
        except requests.HTTPError:
            if i == max_tries - 1: raise
            time.sleep((backoff ** i) + (0.1 * i))
        except requests.RequestException:
            if i == max_tries - 1: raise
            time.sleep((backoff ** i) + (0.1 * i))
    raise RuntimeError("HTTP retry failed")

# ========================= Estruturas =========================
@dataclass
class Region:
    start_1based: int
    end_1based: int
    tag: str
    note: str = ""
    @property
    def length(self) -> int: return self.end_1based - self.start_1based + 1

@dataclass
class PDBUniProtMap:
    uniprot_acc: str
    chain: str
    coverage: float

# ========================= Auditoria BLAST =========================
_BLASTP_AUDIT: List[Dict[str, Any]] = []
_BLASTN_AUDIT: List[Dict[str, Any]] = []
_EXPECTED_TAXID: Optional[int] = None
_EXPECTED_SPECIES: Optional[str] = None

def _species_from_def(s: str) -> str:
    m = re.search(r"\[([^\[\]]+)\]\s*$", s or "")
    return m.group(1) if m else ""

def _audit_hit(kind: str, db: str, region_tag: str, accession: str,
               hit_def: str, pid: float, cov: float, hsp_query_len: int):
    hit_species = _species_from_def(hit_def)
    same_species = None
    if _EXPECTED_SPECIES: same_species = (hit_species == _EXPECTED_SPECIES)
    rec = {
        "kind": kind, "db": db, "region": region_tag, "accession": accession,
        "percent_identity": round(pid * 100.0, 1),
        "query_coverage": round(cov * 100.0, 1),
        "hsp_query_len": int(hsp_query_len),
        "definition": (hit_def or "")[:200],
        "species": hit_species, "same_species": same_species
    }
    ( _BLASTP_AUDIT if kind=="protein" else _BLASTN_AUDIT ).append(rec)
    if CONFIG["output"].get("blast_progress") and CONFIG["output"].get("log_hits", True):
        ss = "" if same_species is None else (" ✓same" if same_species else " ✗off")
        print(f"[HIT {'P' if kind=='protein' else 'N'}] db={db} region={region_tag} acc={accession} "
              f"pid≈{rec['percent_identity']}% cov≈{rec['query_coverage']}% sp='{hit_species}'{ss} "
              f"def={rec['definition']}")

# ========================= UniProt helpers =========================
def fetch_uniprot_json(acc_or_iso: str) -> Dict:
    url = f"https://rest.uniprot.org/uniprotkb/{acc_or_iso}.json"
    return _http_retry("GET", url).json()

def fetch_uniprot_fasta(acc: str, include_isoforms: bool=False) -> Dict[str, str]:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    if include_isoforms: url += "?includeIsoform=true"
    r = _http_retry("GET", url)
    out: Dict[str, str] = {}; cur_id=None; buf=[]
    for ln in r.text.splitlines():
        if ln.startswith(">"):
            if cur_id and buf: out[cur_id] = "".join(buf).strip()
            hdr = ln[1:].strip()
            m = re.search(r"\|([A-Z0-9]+(?:-\d+)?)\|", hdr)
            cur_id = m.group(1) if m else hdr.split()[0]
            buf=[]
        else:
            buf.append(ln.strip())
    if cur_id and buf: out[cur_id] = "".join(buf).strip()
    return out

def uniprot_sequence_from_json(entry: Dict) -> str:
    v = entry.get("sequence", {}).get("value") or ""
    return re.sub(r"\s+", "", v)

def uniprot_taxonomy(entry: Dict) -> Tuple[str, Optional[int]]:
    org = entry.get("organism", {}).get("scientificName") or "?"
    taxid = entry.get("organism", {}).get("taxonId")
    try: taxid = int(taxid) if taxid is not None else None
    except Exception: taxid = None
    return org, taxid

def extract_features_as_regions(entry: Dict, seq_len: int) -> List[Region]:
    cfg = CONFIG["regions"]
    if not cfg.get("use_uniprot_features", True): return []
    incl = set(t.lower() for t in cfg.get("include_feature_types", []))
    point_flank = int(cfg.get("point_flank", 5))
    min_len = int(cfg.get("default_min_len", 1))
    do_merge = bool(cfg.get("merge_overlaps", True))
    out: List[Region] = []
    for feat in entry.get("features", []) or []:
        ftype = (feat.get("type") or "").lower()
        if incl and ftype not in incl: continue
        loc = feat.get("location") or {}
        beg = loc.get("start", {}).get("value"); end = loc.get("end", {}).get("value")
        if beg is None or end is None: continue
        try: beg = int(beg); end = int(end)
        except Exception: continue
        note = (feat.get("description") or "")[:140].strip()
        if beg == end:
            beg = max(1, beg - point_flank); end = min(seq_len, end + point_flank)
        if end >= beg and (end - beg + 1) >= min_len:
            out.append(Region(beg, end, ftype.upper(), note))
    if do_merge and out:
        out.sort(key=lambda r: (r.tag, r.start_1based, r.end_1based))
        merged: List[Region] = []; cur = out[0]
        for r in out[1:]:
            if r.tag == cur.tag and r.start_1based <= cur.end_1based + 1:
                cur.end_1based = max(cur.end_1based, r.end_1based)
            else:
                merged.append(cur); cur = r
        merged.append(cur); out = merged
    if cfg.get("add_full", True):
        out.insert(0, Region(1, seq_len, "FULL", "FULL"))
    return out

# ========================= Gene names / gating =========================
# Mini-mapa de parálogos comuns (evita poluição em EP300/CREBBP)
_PARALOG_EXCLUDES: Dict[str, List[str]] = {
    "EP300":  ["CREBBP","CBP","KAT3A"],
    "CREBBP": ["EP300","P300","KAT3B"],
    "KAT3A":  ["KAT3B","CREBBP","CBP","EP300","P300"],
    "KAT3B":  ["KAT3A","EP300","P300","CREBBP","CBP"],
}

def _collect_gene_names_from_uniprot(entry: Dict) -> Set[str]:
    names: Set[str] = set()
    for g in (entry.get("genes") or []):
        v = (g.get("geneName") or {}).get("value")
        if v: names.add(v)
        for syn in (g.get("synonyms") or []):
            sv = syn.get("value")
            if sv: names.add(sv)
    # Também pega "recommendedName" long name → extrai siglas entre parênteses (ex.: "... (EP300)")
    rec = (entry.get("proteinDescription") or {}).get("recommendedName") or {}
    full = rec.get("fullName", {}).get("value")
    if isinstance(full, str):
        for m in re.findall(r"\(([A-Za-z0-9\-]+)\)", full):
            names.add(m)
    return {n for n in names if re.search(r"[A-Za-z0-9]", n)}

def _collect_gene_names_from_genpept_acc(prot_acc: str) -> Set[str]:
    try:
        handle = Entrez.efetch(db="protein", id=prot_acc, rettype="gp", retmode="text")
        record = SeqIO.read(handle, "genbank"); handle.close()
    except Exception:
        return set()
    names: Set[str] = set()
    for feat in record.features:
        if feat.type in ("CDS","gene"):
            for key in ("gene","gene_synonym"):
                for v in feat.qualifiers.get(key, []):
                    for token in re.split(r"[;,/]\s*", v):
                        if token: names.add(token)
    return {n for n in names if re.search(r"[A-Za-z0-9]", n)}

def _build_gene_gate_patterns(accept_names: Set[str], extra_accept: List[str], exclude_names: List[str]):
    acc = set(accept_names) | set(extra_accept or [])
    exc = set(exclude_names or [])
    # auto-excludes para parálogos comuns
    for k, vals in _PARALOG_EXCLUDES.items():
        if k in acc: exc.update(vals)
    # Regex com borda de palavra; case-insensitive
    accept_re = re.compile(r"\b(" + "|".join(map(re.escape, sorted(acc))) + r")\b", re.I) if acc else None
    exclude_re = re.compile(r"\b(" + "|".join(map(re.escape, sorted(exc))) + r")\b", re.I) if exc else None
    return accept_re, exclude_re

def _hit_title(aln) -> str:
    for attr in ("hit_def","title"):
        v = getattr(aln, attr, None)
        if isinstance(v, str) and v: return v
    return ""

def _gene_gate_allows(title: str, accept_re, exclude_re) -> bool:
    if exclude_re and exclude_re.search(title): return False
    if accept_re and not accept_re.search(title): return False
    return True

# ========================= UniProt ⇄ outras bases =========================
def _strip_refseq_version(acc: str) -> str:
    try: return acc.split('.')[0]
    except Exception: return acc

def map_refseq_protein_to_uniprot(ncbi_prot_acc: str) -> List[str]:
    base = _strip_refseq_version(ncbi_prot_acc)
    tried: List[str] = []
    def q(qs: str) -> List[str]:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": qs, "fields": "accession", "format": "json", "size": "50"}
        tried.append(f"{url}?query={qs}")
        try:
            j = _http_retry("GET", url, params=params).json() or {}
            out = []
            for it in (j.get("results") or []):
                acc = it.get("primaryAccession") or it.get("uniProtkbId")
                if acc: out.append(str(acc))
            return out
        except Exception:
            return []
    results: List[str] = []
    for qs in [f"xref:RefSeq:{ncbi_prot_acc}", f"xref:RefSeq:{base}",
               f'database:RefSeq AND "{ncbi_prot_acc}"', f'database:RefSeq AND "{base}"',
               f"xref:RefSeq_Protein:{ncbi_prot_acc}", f"xref:RefSeq_Protein:{base}",
               f'"{ncbi_prot_acc}"', f'"{base}"', ncbi_prot_acc, base]:
        if results: break
        results += q(qs)
    if not results and CONFIG["output"].get("blast_progress"):
        print(f"[DIAG] UniProt search sem resultados para {ncbi_prot_acc}.")
    return list(dict.fromkeys(results))

# ========================= PDBe SIFTS (PDB → UniProt) =========================
def fetch_pdb_uniprot_mappings(pdb_id: str) -> List[PDBUniProtMap]:
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    j = _http_retry("GET", url).json()
    maps: List[PDBUniProtMap] = []
    block = j.get(pdb_id.lower())
    if not block: return maps
    for acc, obj in (block.get("UniProt") or {}).items():
        seq_len = int(obj.get("sequence_length") or 0) or None
        for m in (obj.get("mappings") or []):
            chain = (m.get("chain_id") or m.get("chain") or "?")
            try: start = int(m.get("unp_start") or 1); end = int(m.get("unp_end") or start)
            except Exception: start, end = 1, 1
            cov = m.get("coverage")
            if cov is None:
                cov = ((end-start+1)/float(seq_len)) if (seq_len and end>=start) else 0.0
            try: cov = max(0.0, min(1.0, float(cov)))
            except Exception: cov = 0.0
            maps.append(PDBUniProtMap(uniprot_acc=acc, chain=str(chain), coverage=cov))
    return maps

def choose_best_mapping(maps: List[PDBUniProtMap]) -> PDBUniProtMap:
    if not maps: raise ValueError("Sem mapeamentos SIFTS")
    return max(maps, key=lambda m: m.coverage)

# ========================= Cross-refs UniProt → nuccore/protein =========================
def extract_nuccore_and_prot_from_uniprot(entry: Dict) -> Tuple[List[str], List[str]]:
    nuccore: List[str] = []; prot_ids: List[str] = []
    for x in entry.get("uniProtKBCrossReferences", []) or []:
        db = (x.get("database") or "").upper(); xid = (x.get("id") or "").strip()
        if db == "REFSEQ" and xid and re.match(r"^[NXYP]P_[0-9]+\.[0-9]+$", xid): prot_ids.append(xid)
        if db in ("REFSEQ","EMBL","GENBANK","DDBJ"):
            for p in x.get("properties", []) or []:
                k = (p.get("key") or "").lower(); v = (p.get("value") or "").strip()
                if not v: continue
                if "nucleotide" in k or "nucleotid" in k:
                    for s in re.split(r"[;,]\s*", v):
                        if re.match(r"^[A-Z]{1,3}[_\.][A-Za-z0-9_\.]+$", s): nuccore.append(s)
    return list(dict.fromkeys(nuccore)), list(dict.fromkeys(prot_ids))

# ========================= GenBank CDS mapping =========================
def elink_protein_to_nuccore_accs(prot_acc: str) -> List[str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {"dbfrom":"protein","db":"nuccore","id":prot_acc,"linkname":"protein_nuccore","retmode":"json"}
    if NCBI_API_KEY: params["api_key"]=NCBI_API_KEY
    if NCBI_EMAIL: params["email"]=NCBI_EMAIL
    r = _http_retry("GET", url, params=params)
    ids: List[str] = []
    try:
        j = r.json()
        for ls in (j.get("linksets") or []):
            for ldb in (ls.get("linksetdbs") or []):
                ids.extend([str(x) for x in (ldb.get("links") or [])])
    except Exception:
        import xml.etree.ElementTree as ET
        root = ET.fromstring(r.content)
        for el in root.findall(".//LinkSetDb/Link/Id"):
            if el is not None and el.text: ids.append(el.text.strip())
    if not ids: return []
    url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    p2 = {"db":"nuccore","id":",".join(ids),"rettype":"acc","retmode":"text"}
    if NCBI_API_KEY: p2["api_key"]=NCBI_API_KEY
    if NCBI_EMAIL: p2["email"]=NCBI_EMAIL
    r2 = _http_retry("GET", url2, params=p2)
    out = [ln.strip() for ln in r2.text.splitlines() if ln.strip()]
    return list(dict.fromkeys(out))

def pick_cds_from_genbank(record, prot_seq: str, candidate_protein_ids: List[str]) -> Tuple[Optional[str], Optional[str], Optional[Dict]]:
    best = None
    for feat in record.features:
        if feat.type != "CDS": continue
        q = feat.qualifiers or {}
        protein_id = (q.get("protein_id",[None])[0]) or None
        translation = (q.get("translation",[None])[0]) or None
        codon_start = int((q.get("codon_start",[1])[0]) or 1)
        transl_table = int((q.get("transl_table",[1])[0]) or 1)

        spliced = str(feat.extract(record.seq)).upper()
        if codon_start in (2,3): spliced = spliced[codon_start-1:]
        if len(spliced) % 3 != 0: spliced = spliced[:len(spliced)//3*3]
        try: tr = str(Seq(spliced).translate(table=transl_table, to_stop=False))
        except Exception: tr = None
        tr_nostop = tr[:-1] if (tr and tr.endswith("*")) else tr

        match_protid = (protein_id in set(candidate_protein_ids)) if protein_id else False
        prot_seq_nostop = prot_seq[:-1] if prot_seq.endswith("*") else prot_seq
        match_translation = (translation == prot_seq) or (tr_nostop == prot_seq_nostop)

        score = (1 if match_protid else 0, 1 if match_translation else 0)
        if score == (0,0): continue
        meta = {"protein_id": protein_id, "codon_start": codon_start, "transl_table": transl_table}
        cand = (spliced, translation, meta, score)
        if best is None or cand[3] > best[3]: best = cand
    if not best: return None, None, None
    return best[0], best[1], best[2]

def find_cds_strict(nuccore_accs: List[str], prot_seq: str, candidate_protein_ids: List[str]) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[Dict]]:
    for acc in nuccore_accs:
        try:
            if CONFIG["output"]["blast_progress"]: print(f"[DIAG] Verificando CDS em {acc} …")
            handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank"); handle.close()
            cds_nt, translation, meta = pick_cds_from_genbank(record, prot_seq, candidate_protein_ids)
            if not cds_nt: continue
            prot_seq_nostop = prot_seq[:-1] if prot_seq.endswith("*") else prot_seq
            if len(cds_nt) != 3*len(prot_seq_nostop):
                tr = str(Seq(cds_nt).translate(to_stop=False))
                if tr.endswith("*"): tr = tr[:-1]
                if tr != prot_seq_nostop: continue
            return acc, cds_nt, translation, meta
        except Exception as e:
            if CONFIG["output"]["blast_progress"]: print(f"[DIAG]  - falha em {acc}: {e}")
            continue
    return None, None, None, None

# ========================= TaxID helpers =========================
def _get_tax_from_genpept(prot_acc: str) -> Tuple[Optional[int], Optional[str]]:
    tax, org = None, None
    try:
        handle = Entrez.efetch(db="protein", id=prot_acc, rettype="gp", retmode="text")
        record = SeqIO.read(handle, "genbank"); handle.close()
    except Exception:
        return None, None
    try: org = record.annotations.get("organism")
    except Exception: org = None
    try:
        for feat in record.features:
            if feat.type == "source":
                for x in feat.qualifiers.get("db_xref", []):
                    if x.startswith("taxon:"): tax = int(x.split(":")[1])
                if not org: org = (feat.qualifiers.get("organism", [""])[0] or None)
    except Exception: pass
    return tax, org

def ensure_blast_same_species_filters_from_taxid(taxid: Optional[int], species: Optional[str]=None):
    global _EXPECTED_TAXID, _EXPECTED_SPECIES
    if species: _EXPECTED_SPECIES = species
    if taxid and CONFIG.get("blast", {}).get("same_species_only", False):
        q = f"txid{int(taxid)}[ORGN]"
        CONFIG["blast"]["protein"]["entrez_query"] = q
        CONFIG["blast"]["nt"]["entrez_query"] = q
        _EXPECTED_TAXID = int(taxid)
        if CONFIG["output"].get("blast_progress"): print(f"[BLAST] Filtro por espécie: {q}")

# ========================= BLAST helpers (gating + união) =========================
def _passes_hsp_thresholds(hsp, L: int, min_cov: float, min_id: float) -> Tuple[bool, float, float, int]:
    qseq = hsp.query; q_aligned = sum(1 for c in qseq if c != '-')
    if q_aligned == 0: return False, 0.0, 0.0, 0
    cov = q_aligned / max(L, 1); 
    if cov < min_cov: return False, cov, 0.0, q_aligned
    ident = hsp.identities / q_aligned
    if ident < min_id: return False, cov, ident, q_aligned
    return True, cov, ident, q_aligned

# -------- BLASTp (janela) --------
def blastp_region(seq_window: str, entrez_query: Optional[str], params: Dict) -> Dict[int, set]:
    if not seq_window or not re.search(r"[A-Z]", seq_window): return {}
    handle = NCBIWWW.qblast(
        program="blastp",
        database=params.get("db", "nr"),
        sequence=seq_window,
        hitlist_size=int(params.get("hitlist_size", 25)),
        expect=float(params.get("expect", 1e-5)),
        entrez_query=entrez_query,
        service="plain",
        megablast=False,
        filter=bool(params.get("filter_low_complexity", True))
    )
    rec = NCBIXML.read(handle); handle.close()
    L = len(seq_window)
    min_id = float(params.get("min_identity", 0.97))
    min_cov = float(params.get("min_query_coverage", 0.95))
    region_tag = params.get("audit_region_tag", "?")
    db_used = params.get("db", "nr")
    accept_re = params.get("gene_accept_re"); exclude_re = params.get("gene_exclude_re")
    relax = bool(params.get("gene_gate_relax", False))

    variants: Dict[int, set] = {}
    accepted_any = False

    def _consume(aln, hsp, audit_db_suffix=""):
        ok, cov, ident, q_aligned = _passes_hsp_thresholds(hsp, L, min_cov, min_id)
        if not ok: return False
        _audit_hit("protein", db_used + audit_db_suffix, region_tag,
                   getattr(aln, "accession", "") or "", _hit_title(aln), ident, cov, q_aligned)
        qpos = 0
        for qc, hc in zip(hsp.query, hsp.sbjct):
            if qc != '-':
                qpos += 1
                if hc != '-' and qc != hc:
                    variants.setdefault(qpos, set()).add(hc)
        return True

    # Passo 1: com gene gating
    for aln in rec.alignments:
        title = _hit_title(aln)
        if not _gene_gate_allows(title, accept_re, exclude_re):
            continue
        for hsp in aln.hsps:
            if _consume(aln, hsp):
                accepted_any = True

    # Passo 2 (opcional): RELAX (ignora gating se nada passou)
    if relax and not accepted_any:
        if CONFIG["output"].get("blast_progress"):
            print(f"[GATE] Sem hits aceitos em {region_tag}. RELAX (sem gene gating).")
        for aln in rec.alignments:
            for hsp in aln.hsps:
                _consume(aln, hsp, audit_db_suffix=" RELAX")

    return variants

def blastp_region_multi(seq_window: str, entrez_query: Optional[str], params: Dict) -> Dict[int, set]:
    dbs = params.get("dbs") or [params.get("db", "nr")]
    union: Dict[int, set] = {}
    for db in dbs:
        if CONFIG["output"]["blast_progress"]:
            print(f"[BLASTp] {params.get('audit_region_tag','?')} em {db}… (len={len(seq_window)} aa)")
        local = blastp_region(seq_window, entrez_query, {**params, "db": db})
        for pos, alts in local.items():
            union.setdefault(pos, set()).update(alts)
    return union

# -------- BLASTn (janela DNA/ RNA) --------
def blastn_region(nt_window: str, db: str, params: Dict, is_mrna: bool=False) -> Dict[int, set]:
    if not nt_window or not re.search(r"[ACGTU]", nt_window, re.I): return {}
    query_dna = nt_window.upper().replace("U", "T")
    entrez_q = params.get("entrez_query")
    if CONFIG["output"]["blast_progress"]:
        print(f"[BLASTn] {params.get('audit_region_tag','?')} em {db} (len={len(query_dna)} nt) "
              f"{'(mesma espécie)' if entrez_q else '(todas as espécies)'}…")
    handle = NCBIWWW.qblast(
        program="blastn",
        database=db,
        sequence=query_dna,
        hitlist_size=int(params.get("hitlist_size", 25)),
        expect=float(params.get("expect", 1e-10)),
        service="plain",
        megablast=bool(params.get("megablast", True)),
        entrez_query=entrez_q,
        filter=bool(params.get("filter_low_complexity", True))
    )
    rec = NCBIXML.read(handle); handle.close()

    L = len(query_dna)
    min_id = float(params.get("min_identity", 0.98))
    min_cov = float(params.get("min_query_coverage", 0.95))
    region_tag = params.get("audit_region_tag", "?")
    accept_re = params.get("gene_accept_re"); exclude_re = params.get("gene_exclude_re")
    relax = bool(params.get("gene_gate_relax", False))

    variants: Dict[int, set] = {}
    accepted_any = False

    def _log_best(rec_obj, db_suffix=""):
        for aln in rec_obj.alignments:
            best_cov = 0.0; best_pid = 0.0; best_hsp = None
            title = _hit_title(aln)
            # log somente dos aceitos (respeitando gating neste nível)
            if not _gene_gate_allows(title, accept_re, exclude_re) and not db_suffix.endswith("RELAX"):
                continue
            for hsp in aln.hsps:
                ok, cov, pid, _qa = _passes_hsp_thresholds(hsp, L, min_cov, min_id)
                if ok and (pid, cov) > (best_pid, best_cov):
                    best_pid, best_cov, best_hsp = pid, cov, hsp
            if best_hsp:
                _audit_hit("nucleotide", f"{db}{db_suffix}", region_tag,
                           getattr(aln, "accession","") or "", title, best_pid, best_cov, best_hsp.align_length)

    def _consume(aln, hsp, db_suffix=""):
        ok, cov, ident, q_aligned = _passes_hsp_thresholds(hsp, L, min_cov, min_id)
        if not ok: return False
        _audit_hit("nucleotide", f"{db}{db_suffix}", region_tag,
                   getattr(aln,"accession","") or "", _hit_title(aln), ident, cov, q_aligned)
        qpos = 0
        for qc, hc in zip(hsp.query.upper(), hsp.sbjct.upper()):
            if qc != '-':
                qpos += 1
                if hc != '-' and qc != hc:
                    variants.setdefault(qpos, set()).add(hc)
        return True

    # 1) com gene gating
    if CONFIG["output"]["blast_progress"]: _log_best(rec, "")
    for aln in rec.alignments:
        if not _gene_gate_allows(_hit_title(aln), accept_re, exclude_re):
            continue
        for hsp in aln.hsps:
            if _consume(aln, hsp):
                accepted_any = True

    # 2) RELAX, se necessário
    if relax and not accepted_any:
        if CONFIG["output"]["blast_progress"]:
            print(f"[GATE] Sem hits aceitos em NT {region_tag}. RELAX (sem gene gating).")
        if CONFIG["output"]["blast_progress"]: _log_best(rec, " RELAX")
        for aln in rec.alignments:
            for hsp in aln.hsps:
                _consume(aln, hsp, db_suffix=" RELAX")

    # 3) Preview off-species (somente log)
    try:
        if CONFIG["blast"].get("log_offspecies_preview", False) and entrez_q:
            if CONFIG["output"]["blast_progress"]: print("[BLASTn] Preview off-species (log)…")
            handle2 = NCBIWWW.qblast(
                program="blastn", database=db, sequence=query_dna,
                hitlist_size=min(int(params.get("hitlist_size", 25)), 5),
                expect=float(params.get("expect", 1e-10)), service="plain",
                megablast=bool(params.get("megablast", True)), entrez_query=None,
                filter=bool(params.get("filter_low_complexity", True))
            )
            rec2 = NCBIXML.read(handle2); handle2.close()
            _log_best(rec2, " PREVIEW")
    except Exception as e:
        if CONFIG["output"]["blast_progress"]:
            print(f"[WARN] Preview off-species falhou: {e}")

    return variants

# -------- BLASTp FULL (agrega variantes absolutas) --------
def _aggregate_variants_from_full_hsp(seq_len: int, hsp, aln, db_used: str, region_label: str,
                                      min_id: float, min_cov: float,
                                      variants_abs: Dict[int, set], covered_abs: set):
    ok, cov, ident, q_aligned = _passes_hsp_thresholds(hsp, seq_len, min_cov, min_id)
    if not ok: return
    _audit_hit("protein", db_used, region_label,
               getattr(aln,"accession","") or "", _hit_title(aln), ident, cov, q_aligned)
    qabs = int(getattr(hsp, "query_start", 1)) - 1
    for qc, hc in zip(hsp.query, hsp.sbjct):
        if qc != '-':
            qabs += 1
            if 1 <= qabs <= seq_len:
                covered_abs.add(qabs)
                if hc != '-' and qc != hc:
                    variants_abs.setdefault(qabs, set()).add(hc)

def blastp_full_aggregate_variants(seq_full: str, entrez_query: Optional[str], params: Dict) -> Tuple[Dict[int, set], set]:
    if not seq_full or not re.search(r"[A-Z]", seq_full): return {}, set()
    L = len(seq_full)
    variants_abs: Dict[int, set] = {}; covered_abs: set = set()
    dbs = params.get("full_dbs") or params.get("dbs") or [params.get("db", "nr")]
    min_id = float(params.get("min_identity", 0.85))
    min_cov_full = float(params.get("full_min_query_coverage", params.get("min_query_coverage", 0.90)))
    hitlist_full = int(params.get("full_hitlist_size", min(params.get("hitlist_size", 200), 50)))
    max_align = int(params.get("max_alignments_process", 200))
    accept_re = params.get("gene_accept_re"); exclude_re = params.get("gene_exclude_re")
    relax = bool(params.get("gene_gate_relax", False))
    accepted_any = False

    for db in dbs:
        if CONFIG["output"]["blast_progress"]:
            print(f"[BLASTp] FULL ({L} aa) em {db} (entrez_query={params.get('entrez_query')})…")
        handle = NCBIWWW.qblast(
            program="blastp", database=db, sequence=seq_full,
            hitlist_size=hitlist_full, expect=float(params.get("expect", 1e-5)),
            entrez_query=params.get("entrez_query"), service="plain", megablast=False,
            filter=bool(params.get("filter_low_complexity", True))
        )
        rec = NCBIXML.read(handle); handle.close()

        # 1) com gene gating
        count_aln = 0
        for aln in rec.alignments:
            count_aln += 1
            if count_aln > max_align: break
            if not _gene_gate_allows(_hit_title(aln), accept_re, exclude_re): continue
            for hsp in aln.hsps:
                _aggregate_variants_from_full_hsp(L, hsp, aln, db, "FULL",
                                                  min_id, min_cov_full, variants_abs, covered_abs)
                accepted_any = True

        # 2) RELAX se nada aceito
        if relax and not accepted_any:
            if CONFIG["output"]["blast_progress"]:
                print("[GATE] FULL sem hits aceitos. RELAX (sem gene gating).")
            count_aln = 0
            for aln in rec.alignments:
                count_aln += 1
                if count_aln > max_align: break
                for hsp in aln.hsps:
                    _aggregate_variants_from_full_hsp(L, hsp, aln, db+" RELAX", "FULL",
                                                      min_id, min_cov_full, variants_abs, covered_abs)

    return variants_abs, covered_abs

# ========================= Effatha builders =========================
def effatha_aa(seq: str, region: Region, local_variants: Dict[int, set]) -> str:
    s = region.start_1based; e = region.end_1based; buf=[]
    for i, pos in enumerate(range(s, e+1), start=1):
        ref = seq[pos-1]; alts = sorted(x for x in local_variants.get(i, set()) if x != ref)
        buf.append("[" + "/".join([ref] + alts) + "]" if alts else ref)
    return "".join(buf)

def effatha_nt(nt_seq: str, nt_variants: Dict[int, set], use_u: bool=False) -> str:
    buf=[]
    for i, ref in enumerate(nt_seq, start=1):
        rc = 'U' if (use_u and ref.upper()=='T') else ref
        alts = sorted(('U' if (use_u and a.upper()=='T') else a) for a in nt_variants.get(i, set()) if a != rc)
        buf.append("[" + "/".join([rc] + alts) + "]" if alts else rc)
    return "".join(buf)

# ========================= Núcleo =========================
def core_pipeline_using_uniprot(target_uniprot: str,
                                refseq_hint_seq: Optional[str]=None,
                                mapping_info: Optional[Dict]=None,
                                extra_prot_ids: Optional[List[str]]=None):
    print(f"[INFO] UniProt {target_uniprot}: coletando sequência, features e CDS…")
    entry = fetch_uniprot_json(target_uniprot)
    seq = uniprot_sequence_from_json(entry)
    org, taxid = uniprot_taxonomy(entry)

    ensure_blast_same_species_filters_from_taxid(taxid, species=org)

    # ---- Gene gating patterns
    gate_cfg = CONFIG["blast"]["gene_gate"]
    accept_names = _collect_gene_names_from_uniprot(entry) if gate_cfg.get("accept_from_uniprot", True) else set()
    accept_re, exclude_re = _build_gene_gate_patterns(accept_names, gate_cfg.get("extra_accept", []), gate_cfg.get("exclude", []))

    regs = extract_features_as_regions(entry, len(seq))

    if CONFIG["regions"].get("fallback_genpept_if_uniprot_featureless", True):
        if (not regs) or all(r.tag == "FULL" for r in regs):
            _nuccs_from_up, prot_ids_from_up = extract_nuccore_and_prot_from_uniprot(entry)
            prot_candidates = list(dict.fromkeys((extra_prot_ids or []) + prot_ids_from_up))
            if prot_candidates:
                print(f"[INFO] UniProt sem features. Fallback GenPept usando {prot_candidates[0]} …")
                try:
                    regs_gp = extract_regions_from_genpept_features(prot_candidates[0])
                    if regs_gp:
                        regs = regs_gp
                        print("[OK] Features do GenPept aplicadas.")
                except Exception as e:
                    print(f"[WARN] Fallback GenPept falhou: {e}")
    if not regs: regs = [Region(1, len(seq), "FULL", "FULL")]

    # ===== AA (FULL + regiões) =====
    aa_by_region: List[Tuple[str,int,int,str,str]] = []
    full_variants_abs: Dict[int, set] = {}; full_covered_abs: set = set()
    if CONFIG["blast"]["enable"]:
        bp = {**CONFIG["blast"]["protein"],
              "gene_accept_re": accept_re, "gene_exclude_re": exclude_re,
              "gene_gate_relax": bool(CONFIG["blast"]["gene_gate"]["relax_if_no_hits"])}
        entrez_q = bp.get("entrez_query")

        if bp.get("smart_full_first", True):
            full_variants_abs, full_covered_abs = blastp_full_aggregate_variants(seq, entrez_q, bp)

        for r in regs:
            frag = seq[r.start_1based-1:r.end_1based]
            local_vars: Dict[int, set] = {}
            if full_variants_abs:
                for i_abs in range(r.start_1based, r.end_1based + 1):
                    if i_abs in full_variants_abs:
                        local_vars.setdefault(i_abs - r.start_1based + 1, set()).update(full_variants_abs[i_abs])

            need_window_blast = False
            if bp.get("smart_full_first", True):
                region_len = r.length
                if region_len >= int(bp.get("regions_min_len_blast", 40)):
                    region_positions = set(range(r.start_1based, r.end_1based + 1))
                    covered_here = region_positions & full_covered_abs
                    gap_frac = 1.0 - (len(covered_here) / float(region_len)) if region_len else 1.0
                    skip_tags = set(t.strip().upper() for t in bp.get("skip_tags_for_region_blast", []))
                    need_window_blast = (gap_frac > float(bp.get("gap_frac_threshold", 0.15))) and ((r.tag or "").upper() not in skip_tags)
                    if CONFIG["output"]["blast_progress"]:
                        print(f"[BLASTp] {r.tag} ({region_len} aa) — cobertura FULL≈{(1-gap_frac):.1%} | fallback_janela={'SIM' if need_window_blast else 'NÃO'}")

            if need_window_blast or not bp.get("smart_full_first", True):
                if CONFIG["output"]["blast_progress"]:
                    print(f"[BLASTp] Fallback janela {r.tag} (entrez_query={entrez_q})…")
                try:
                    bp_with_tag = {**bp, "audit_region_tag": r.tag}
                    local_extra = blastp_region_multi(frag, entrez_q, bp_with_tag)
                    for k, alts in local_extra.items():
                        local_vars.setdefault(k, set()).update(alts)
                except Exception as e:
                    if CONFIG["output"]["blast_progress"]:
                        print(f"[WARN] BLASTP falhou em {r.tag} {r.start_1based}-{r.end_1based}: {e}")
            aa_by_region.append((r.tag, r.start_1based, r.end_1based, effatha_aa(seq, r, local_vars), r.note))
    else:
        for r in regs:
            aa_by_region.append((r.tag, r.start_1based, r.end_1based, seq[r.start_1based-1:r.end_1based], r.note))

    # ===== NT por região (CDS real + BLASTn) =====
    nt_segments: List[Tuple[str,int,int,str,str,str]] = []
    if CONFIG["cds_mapping"]["enable"]:
        nuccore_cands, prot_ids_from_up = extract_nuccore_and_prot_from_uniprot(entry)
        extra_from_prot = []
        for pid in (extra_prot_ids or []): extra_from_prot += elink_protein_to_nuccore_accs(pid)
        for pid in (prot_ids_from_up or []): extra_from_prot += elink_protein_to_nuccore_accs(pid)
        nuccore_all = list(dict.fromkeys(nuccore_cands + extra_from_prot))
        if CONFIG["output"]["blast_progress"]:
            print(f"[DIAG] nuccore candidatos: {nuccore_all or '—'} ; protein_ids: {(extra_prot_ids or []) + prot_ids_from_up or ['—']}")

        picked_nuccore, cds_nt, _translation, _meta = find_cds_strict(nuccore_all, seq, (extra_prot_ids or []) + prot_ids_from_up)
        if picked_nuccore and cds_nt:
            bn = {**CONFIG["blast"]["nt"],
                  "gene_accept_re": accept_re, "gene_exclude_re": exclude_re,
                  "gene_gate_relax": bool(CONFIG["blast"]["gene_gate"]["relax_if_no_hits"])}
            for r in regs:
                dna = cds_nt[(r.start_1based-1)*3 : r.end_1based*3]
                mrna = dna.replace("T","U")
                dna_vars = {}; mrna_vars = {}
                if CONFIG["blast"]["enable"]:
                    try:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[BLASTn] DNA {r.tag} ({len(dna)} nt) …")
                        bn_with_tag = {**bn, "audit_region_tag": r.tag}
                        dna_vars = blastn_region(dna, bn.get("dna_db","nt"), bn_with_tag, is_mrna=False)
                    except Exception as e:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[WARN] BLASTN DNA falhou em {r.tag}: {e}")
                    try:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[BLASTn] mRNA {r.tag} ({len(mrna)} nt) …")
                        bn_with_tag = {**bn, "audit_region_tag": r.tag}
                        mrna_vars = blastn_region(mrna, bn.get("rna_db","refseq_rna"), bn_with_tag, is_mrna=True)
                    except Exception as e:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[WARN] BLASTN mRNA falhou em {r.tag}: {e}")
                nt_segments.append((r.tag, r.start_1based, r.end_1based,
                                    effatha_nt(dna, dna_vars, use_u=False),
                                    effatha_nt(mrna, mrna_vars, use_u=True),
                                    r.note))
        else:
            if CONFIG["output"]["blast_progress"]:
                print("[DIAG] CDS real não mapeada — NT ausente.")

    protein_meta = {"accession": target_uniprot, "length": len(seq), "organism": org, "taxid": taxid}
    if mapping_info: protein_meta["mapping"] = mapping_info

    region_cards = [{
        "tag": r.tag, "start": r.start_1based, "end": r.end_1based, "length": r.length, "note": r.note,
        "ligand_name": None, "ligand_chebi": None, "ligand_is_metal": None,
        "pct_var_obs": 0.0, "pct_var_fil": 0.0, "max_af": None, "clin_counts": {},
        "pos_with_manual_ev": 0, "sift_tol": 0, "sift_del": 0, "poly_benign": 0, "poly_damaging": 0,
        "has_deletion": False,
    } for r in regs]

    consolidated = {
        "regions": [{"tag": r.tag, "start": r.start_1based, "end": r.end_1based} for r in regs],
        "region_cards": region_cards,
        "aa_regions": [{"tag": t, "start": s, "end": e, "aa": a} for (t,s,e,a,_) in aa_by_region],
        "aa_regions_obs": [], "aa_regions_filtered": [],
        "nt_segments": [{"tag": t, "start": s, "end": e, "dna": dna, "mrna": mrna, "note": note}
                        for (t,s,e,dna,mrna,note) in nt_segments],
    }
    return seq, protein_meta, consolidated

# ========================= Entradas =========================
def run_from_uniprot(uniprot_acc: str):
    seq, meta, consolidated = core_pipeline_using_uniprot(uniprot_acc)
    context = {
        "protein": meta, "layers": {
            "uniprot_variant_enabled": True, "proteins_variation_enabled": True, "blast_enabled": bool(CONFIG["blast"]["enable"])
        },
        "blast": {"variant_positions_blast": 0, "variant_positions_uniprot": 0, "variant_positions_proteins": 0, "variant_positions_filtered": 0},
        **consolidated,
        "sources": {
            "uniprot": consolidated, "pdb": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_protein": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": consolidated["nt_segments"]},
        }
    }
    _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)

def run_from_pdb(pdb_id: str):
    pdb = (pdb_id or "").strip()
    if not pdb: print("[ERRO] PDB vazio."); return
    print(f"Descobrindo UniProt para PDB {pdb} via SIFTS…")
    maps = fetch_pdb_uniprot_mappings(pdb)
    if not maps: print("[ERRO] Nenhum mapeamento via SIFTS."); return
    best = choose_best_mapping(maps); chosen = best.uniprot_acc
    print(f"Selecionado UniProt {chosen} (cadeia {best.chain}, cobertura≈{best.coverage:.1%})")
    org_tax_entry = fetch_uniprot_json(chosen); _org, _taxid = uniprot_taxonomy(org_tax_entry)
    ensure_blast_same_species_filters_from_taxid(_taxid, species=_org)
    mapping_info = {"input_source":"pdb","input_id":pdb_id,"mapped_uniprot":chosen,"notes":"SIFTS (PDBe)"}
    seq, meta, consolidated = core_pipeline_using_uniprot(chosen, mapping_info=mapping_info)
    context = {
        "protein": meta, "layers": {"uniprot_variant_enabled": True, "proteins_variation_enabled": True, "blast_enabled": bool(CONFIG["blast"]["enable"])},
        "blast": {"variant_positions_blast": 0, "variant_positions_uniprot": 0, "variant_positions_proteins": 0, "variant_positions_filtered": 0},
        **consolidated,
        "sources": {
            "uniprot": consolidated, "pdb": consolidated,
            "ncbi_protein": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": consolidated["nt_segments"]},
        }
    }
    _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)

def run_from_ncbi_protein(refseq_acc: str):
    acc = (refseq_acc or "").strip()
    if not acc: raise ValueError("RefSeq Protein ACC vazio.")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db":"protein","id":acc,"rettype":"fasta","retmode":"text"}
    if NCBI_API_KEY: params["api_key"]=NCBI_API_KEY
    if NCBI_EMAIL: params["email"]=NCBI_EMAIL
    r = _http_retry("GET", url, params=params)
    refseq_seq = "".join([ln.strip() for ln in r.text.splitlines() if not ln.startswith(">")])
    if not refseq_seq: raise RuntimeError("Falha ao obter sequência proteica do RefSeq.")

    gp_tax, gp_org = _get_tax_from_genpept(acc)
    ensure_blast_same_species_filters_from_taxid(gp_tax, species=gp_org)

    up_accs = map_refseq_protein_to_uniprot(acc)
    if up_accs:
        # Vai pelo caminho UniProt (com gene gating do UniProt)
        up = choose_best_uniprot_isoform(up_accs, isoform_policy="longest", refseq_hint_seq=refseq_seq)
        seq, meta, consolidated = core_pipeline_using_uniprot(
            up,
            refseq_hint_seq=refseq_seq,
            mapping_info={"input_source":"ncbi_protein","input_id":acc,"mapped_uniprot":up},
            extra_prot_ids=[acc]
        )
        context = {
            "protein": meta, "layers": {"uniprot_variant_enabled": True, "proteins_variation_enabled": True, "blast_enabled": bool(CONFIG["blast"]["enable"])},
            "blast": {"variant_positions_blast": 0, "variant_positions_uniprot": 0, "variant_positions_proteins": 0, "variant_positions_filtered": 0},
            **consolidated,
            "sources": {
                "uniprot": consolidated,
                "pdb": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
                "ncbi_protein": consolidated,
                "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": consolidated["nt_segments"]},
            }
        }
        _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)
        return

    # Fallback: sem UniProt — gene gating a partir do GenPept
    print("[AVISO] Não mapeou RefSeq→UniProt. Usando features + nomes de gene do GenPept.")
    regs = extract_regions_from_genpept_features(acc)
    seq = refseq_seq
    accept_names_gp = _collect_gene_names_from_genpept_acc(acc)
    gate_cfg = CONFIG["blast"]["gene_gate"]
    accept_re, exclude_re = _build_gene_gate_patterns(accept_names_gp, gate_cfg.get("extra_accept", []), gate_cfg.get("exclude", []))

    aa_by_region = []
    if CONFIG["blast"]["enable"]:
        bp = {**CONFIG["blast"]["protein"],
              "gene_accept_re": accept_re, "gene_exclude_re": exclude_re,
              "gene_gate_relax": bool(CONFIG["blast"]["gene_gate"]["relax_if_no_hits"])}
        full_vars_abs, full_cov_abs = ({}, set())
        if bp.get("smart_full_first", True):
            full_vars_abs, full_cov_abs = blastp_full_aggregate_variants(seq, bp.get("entrez_query"), bp)

        for rgn in regs:
            frag = seq[rgn.start_1based-1:rgn.end_1based]; local_vars={}
            if full_vars_abs:
                for i_abs in range(rgn.start_1based, rgn.end_1based + 1):
                    if i_abs in full_vars_abs:
                        local_vars.setdefault(i_abs - rgn.start_1based + 1, set()).update(full_vars_abs[i_abs])
            need_window = True
            if bp.get("smart_full_first", True):
                region_positions = set(range(rgn.start_1based, rgn.end_1based + 1))
                gap = 1.0 - (len(region_positions & full_cov_abs) / float(rgn.length)) if rgn.length else 1.0
                need_window = (rgn.length >= int(bp.get("regions_min_len_blast", 40))) and (gap > float(bp.get("gap_frac_threshold", 0.15))) and ((rgn.tag or "").upper() not in set(t.upper() for t in bp.get("skip_tags_for_region_blast", [])))
            if need_window or not bp.get("smart_full_first", True):
                if CONFIG["output"]["blast_progress"]:
                    print(f"[BLASTp] {rgn.tag} ({rgn.length} aa) …")
                try:
                    bp_with_tag = {**bp, "audit_region_tag": rgn.tag}
                    local_extra = blastp_region_multi(frag, bp.get("entrez_query"), bp_with_tag)
                    for k, alts in local_extra.items():
                        local_vars.setdefault(k, set()).update(alts)
                except Exception as e:
                    if CONFIG["output"]["blast_progress"]:
                        print(f"[WARN] BLASTP falhou: {e}")
            aa_by_region.append({"tag":rgn.tag,"start":rgn.start_1based,"end":rgn.end_1based,"aa": effatha_aa(seq, rgn, local_vars)})
    else:
        for rgn in regs:
            aa_by_region.append({"tag":rgn.tag,"start":rgn.start_1based,"end":rgn.end_1based,"aa": seq[rgn.start_1based-1:rgn.end_1based]})

    nt_segments = []
    if CONFIG["cds_mapping"]["enable"]:
        picked, cds_nt, _t, _m = find_cds_strict(elink_protein_to_nuccore_accs(acc), seq, [acc])
        if picked and cds_nt:
            bn = {**CONFIG["blast"]["nt"],
                  "gene_accept_re": accept_re, "gene_exclude_re": exclude_re,
                  "gene_gate_relax": bool(CONFIG["blast"]["gene_gate"]["relax_if_no_hits"])}
            for rgn in regs:
                dna = cds_nt[(rgn.start_1based-1)*3 : rgn.end_1based*3]; mrna = dna.replace("T","U")
                try:
                    if CONFIG["blast"]["enable"] and CONFIG["output"]["blast_progress"]:
                        print(f"[BLASTn] DNA {rgn.tag} ({len(dna)} nt) …")
                    bn_with_tag = {**bn, "audit_region_tag": rgn.tag}
                    dna_vars = blastn_region(dna, bn.get("dna_db","nt"), bn_with_tag, is_mrna=False) if CONFIG["blast"]["enable"] else {}
                except Exception as e:
                    dna_vars = {}; print(f"[WARN] BLASTN DNA falhou: {e}")
                try:
                    if CONFIG["blast"]["enable"] and CONFIG["output"]["blast_progress"]:
                        print(f"[BLASTn] mRNA {rgn.tag} ({len(mrna)} nt) …")
                    bn_with_tag = {**bn, "audit_region_tag": rgn.tag}
                    mrna_vars = blastn_region(mrna, bn.get("rna_db","refseq_rna"), bn_with_tag, is_mrna=True) if CONFIG["blast"]["enable"] else {}
                except Exception as e:
                    mrna_vars = {}; print(f"[WARN] BLASTN mRNA falhou: {e}")
                nt_segments.append({"tag":rgn.tag,"start":rgn.start_1based,"end":rgn.end_1based,
                                    "dna": effatha_nt(dna, dna_vars, use_u=False),
                                    "mrna": effatha_nt(mrna, mrna_vars, use_u=True),
                                    "note": rgn.note})

    meta = {"accession": acc, "length": len(seq), "organism": _EXPECTED_SPECIES or "?", "taxid": _EXPECTED_TAXID}
    consolidated = {
        "regions": [{"tag":r.tag,"start":r.start_1based,"end":r.end_1based} for r in regs],
        "region_cards": [{
            "tag":r.tag,"start":r.start_1based,"end":r.end_1based,"length":r.length,"note":r.note,
            "ligand_name":None,"ligand_chebi":None,"ligand_is_metal":None,
            "pct_var_obs":0.0,"pct_var_fil":0.0,"max_af":None,"clin_counts":{},
            "pos_with_manual_ev":0,"sift_tol":0,"sift_del":0,"poly_benign":0,"poly_damaging":0,"has_deletion":False
        } for r in regs],
        "aa_regions": aa_by_region, "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": nt_segments,
    }
    context = {
        "protein": meta,
        "layers": {"uniprot_variant_enabled": True, "proteins_variation_enabled": True, "blast_enabled": bool(CONFIG["blast"]["enable"])},
        "blast": {"variant_positions_blast": 0, "variant_positions_uniprot": 0, "variant_positions_proteins": 0, "variant_positions_filtered": 0},
        **consolidated,
        "sources": {
            "uniprot": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "pdb": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_protein": consolidated,
            "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": nt_segments},
        }
    }
    _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)

def run_from_ncbi_gene(gene_cfg: Dict[str, Any]):
    id_type = (gene_cfg or {}).get("id_type") or "entrez"
    taxid = int((gene_cfg or {}).get("taxid") or 9606)
    gene_id = (gene_cfg or {}).get("id") or ""
    symbol = (gene_cfg or {}).get("symbol") or ""
    if id_type == "entrez" and gene_id:
        q = f"{gene_id}[Gene ID] AND txid{taxid}[Organism]"
    elif id_type == "symbol" and symbol:
        q = f"{symbol}[Gene Name] AND txid{taxid}[Organism]"
    else:
        raise ValueError("Config de gene inválida.")

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    p = {"db":"protein","term":q,"retmax":"100","sort":"slen","retmode":"json"}
    if NCBI_API_KEY: p["api_key"]=NCBI_API_KEY
    if NCBI_EMAIL:   p["email"]=NCBI_EMAIL
    r = _http_retry("GET", url, params=p)
    j = r.json()
    ids = (j.get("esearchresult",{}).get("idlist") or [])
    if not ids: raise RuntimeError("Nenhum RefSeq Protein encontrado.")

    url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    p2 = {"db":"protein","id":",".join(ids),"rettype":"acc","retmode":"text"}
    if NCBI_API_KEY: p2["api_key"]=NCBI_API_KEY
    if NCBI_EMAIL:   p2["email"]=NCBI_EMAIL
    r2 = _http_retry("GET", url2, params=p2)
    accs = [ln.strip() for ln in r2.text.splitlines() if ln.strip()]
    if not accs: raise RuntimeError("Falha ao obter ACCs de protein.")
    refseq_acc = max(accs, key=lambda a: len(a))
    return run_from_ncbi_protein(refseq_acc)

# ========================= Saídas =========================
def _write_artifacts_and_context(prot_seq: str, aa_regions: List[Dict[str, Any]], nt_segments: List[Dict[str, Any]], context: Dict[str, Any]):
    rp = CONFIG["output"]["report_txt"]; os.makedirs(os.path.dirname(rp), exist_ok=True)
    with io.open(rp, "w", encoding="utf-8") as fh:
        fh.write("# Effatha — Gene & Protein Analyzer (regiões funcionais)\n\n")
        pmeta = context.get("protein", {})
        fh.write(f"Proteína: {pmeta.get('accession','?')}\n")
        fh.write(f"Organismo: {pmeta.get('organism','?')} (taxid={pmeta.get('taxid','?')})\n")
        fh.write(f"Tamanho (AA): {len(prot_seq)} aa\n\n")
        fh.write("## Aminoácidos por região (sintaxe Effatha)\n")
        for a in aa_regions:
            fh.write(f"\n### {a['tag']} {a['start']}-{a['end']}\n{a['aa']}\n")
        fh.write("\n\n## Segmentos NT (DNA/mRNA; sintaxe Effatha por base)\n")
        if nt_segments:
            for s in nt_segments:
                fh.write(f"\n### {s['tag']} {s['start']}-{s['end']}\nDNA:\n{s['dna']}\n")
                fh.write(f"mRNA:\n{s['mrna']}\n")
        else:
            fh.write("(Sem NT — CDS real não mapeado.)\n")

    rcsv = CONFIG["output"]["export_regions_csv"]
    if rcsv:
        try:
            os.makedirs(os.path.dirname(rcsv), exist_ok=True)
            with io.open(rcsv, "w", encoding="utf-8", newline="") as fcsv:
                w = csv.writer(fcsv); w.writerow(["Tag","Start","End","Length"])
                for a in aa_regions: w.writerow([a["tag"], a["start"], a["end"], a["end"]-a["start"]+1])
        except Exception as e:
            print(f"[WARN] Falha ao escrever regions.csv: {e}")

    vcsv = CONFIG["output"]["export_csv"]
    if vcsv and os.path.basename(vcsv):
        try:
            os.makedirs(os.path.dirname(vcsv), exist_ok=True)
            with io.open(vcsv, "w", encoding="utf-8", newline="") as fcsv:
                w = csv.writer(fcsv); w.writerow(["level","region_tag","pos","observed_set","source"])
        except Exception as e:
            print(f"[WARN] Falha ao escrever variants_blast.csv: {e}")

    hits_csv = CONFIG["output"].get("blast_hits_csv")
    if hits_csv:
        try:
            os.makedirs(os.path.dirname(hits_csv), exist_ok=True)
            with io.open(hits_csv, "w", encoding="utf-8", newline="") as fcsv:
                w = csv.writer(fcsv)
                w.writerow(["kind","db","region","accession","percent_identity","query_coverage","hsp_query_len","species","same_species","definition"])
                for rec in (_BLASTP_AUDIT + _BLASTN_AUDIT):
                    w.writerow([rec.get("kind"),rec.get("db"),rec.get("region"),rec.get("accession"),
                                rec.get("percent_identity"),rec.get("query_coverage"),rec.get("hsp_query_len"),
                                rec.get("species",""),rec.get("same_species",""),rec.get("definition","")])
            if CONFIG["output"].get("blast_progress"): print(f"[OK] BLAST hits salvos em {hits_csv}")
        except Exception as e:
            print(f"[WARN] Falha ao escrever blast_hits_csv: {e}")

    ctx_path = os.path.join(CONFIG["output"]["artifacts_dir"], "context_summary.json")
    try:
        os.makedirs(CONFIG["output"]["artifacts_dir"], exist_ok=True)
        with io.open(ctx_path, "w", encoding="utf-8") as fh:
            json.dump(context, fh, ensure_ascii=False, indent=2)
        print(f"[OK] context_summary.json salvo em {ctx_path}")
    except Exception as e:
        print(f"[ERRO] Ao salvar context_summary.json: {e}")

# ========================= Entrypoint =========================
if __name__ == "__main__":
    _install_signal_handlers(); _start_watchdog()
    try:
        if not NCBI_EMAIL: print("[AVISO] Defina NCBI_EMAIL para cumprir a política do NCBI.")
        src = (CONFIG["input"]["source"] or "").lower()
        if src == "uniprot":
            run_from_uniprot(CONFIG["input"]["uniprot_acc"])
        elif src == "pdb":
            run_from_pdb(CONFIG["input"]["pdb_id"])
        elif src == "ncbi_protein":
            run_from_ncbi_protein(CONFIG["input"]["ncbi_protein_acc"])
        elif src == "ncbi_gene":
            run_from_ncbi_gene(CONFIG["input"]["gene"])
        else:
            print("Defina CONFIG['input']['source'] para 'uniprot' | 'pdb' | 'ncbi_protein' | 'ncbi_gene'.")
        print("[shim] finished main", flush=True)
        _cancel_watchdog(); _finalize_and_exit(0)
    except requests.HTTPError as e:
        _cancel_watchdog(); print(f"[ERRO HTTP] {e} — URL: {getattr(e.response, 'url', '?')}"); _finalize_and_exit(1)
    except SystemExit:
        _cancel_watchdog(); raise
    except Exception as e:
        _cancel_watchdog(); print(f"[ERRO] {e}"); _finalize_and_exit(1)

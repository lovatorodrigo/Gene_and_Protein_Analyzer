from __future__ import annotations
from io import StringIO
from Bio import SeqIO

# ==============================
# CONFIG — ÚNICO PONTO DE AJUSTE
# ==============================
CONFIG = {
    "input": {
        # "pdb" | "uniprot" | "ncbi_protein" | "ncbi_gene"
        "source": "ncbi_protein",
        "pdb_id": "5xnl",
        "uniprot_acc": "O14744",
        "ncbi_protein_acc": "AML61188.1",
        "gene": {
            "id_type": "entrez",      # "entrez" | "symbol"
            "id": "7157",             # GeneID se id_type=="entrez"
            "symbol": "TP53",         # símbolo se id_type=="symbol"
            "taxid": 9606,            # obrigatório p/ symbol
            "isoform_policy": "longest"  # "longest" | "mane" | "all"
        },
    },

    "regions": {
        "mode": "uniprot_features",      # "uniprot_features" | "pfam" | "both"
        "feature_types": "ALL",
        "flank_left": 5,
        "flank_right": 5,
        "merge_overlaps": True,
        # incluir features NCBI (GenPept) quando houver RefSeq associado
        "include_ncbi_protein_features": True,
    },

    "windows": {
        "flank_scan": 5,
        "stride": 10
    },

    "blast": {
        "enable": True,
        "database": "nr",
        "expect": 1e-5,
        "hitlist_size": 10,
        "throttle_seconds": 11,
        "min_pct_identity": 90.0,
        "min_query_coverage": 0.95,
        "same_species_only": True
    },

    "output": {
        "blast_progress": True,
        "export_csv": None,
        "report_txt": "run_report.txt",
        "export_regions_csv": None,
        "artifacts_dir": ".",
    },

    "cds_mapping": {
        "enable": True
    },

    "external_layers": {
        "uniprot_variant": {"enable": True},
        "variation_api": {"enable": True},
    },

    "variation_filters": {
        "enabled": True,
        "min_af": 0.01,
        "clinical": {
            "allow": ["benign","likely_benign"],
            "allow_if_conflicting": False,
            "allow_uncertain": False,
        },
        "evidence": {
            "assertion": "prefer_manual",
            "treat_missing_as": "include",
        },
        "predictions": {
            "sift_allow": ["tolerated","tolerated_low_confidence"],
            "polyphen_allow": ["benign","possibly_benign"],
            "require_both": False,
        },
        "include_blast_in_filtered": True,
    },

    "net": {
        "user_agent": "Effatha-RnD/1.0 (contact: lovato.rodrigo@gmail.com)",
        "timeout_connect": 10,
        "timeout_read": 30,
        "ncbi_api_key_env": "NCBI_API_KEY",
        "eutils_email_env": "EUTILS_EMAIL"
    },
}
# ==============================
# FIM DO BLOCO CONFIG
# ==============================

from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Iterable, Set
import argparse
import os
import time
import logging
import csv
import json
from datetime import datetime

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from io import StringIO
import xml.etree.ElementTree as ET

# ----------------------------------------------------------------------------
# Lista completa de feature types da UniProt (atual)
# ----------------------------------------------------------------------------
ALL_UNIPROT_FEATURE_TYPES = [
    "ACT_SITE","BINDING","CARBOHYD","CHAIN","COILED","COMPBIAS","CONFLICT","CROSSLNK",
    "DNA_BIND","DISULFID","DOMAIN","HELIX","INIT_MET","INTRAMEM","LIPID","MOD_RES",
    "MOTIF","MUTAGEN","NON_CONS","NON_STD","NON_TER","PEPTIDE","PROPEP","REGION",
    "REPEAT","SIGNAL","SITE","STRAND","TOPO_DOM","TRANSIT","TRANSMEM","TURN",
    "VARIANT","VAR_SEQ",
]

def _resolve_feature_types(cfg_value):
    if isinstance(cfg_value, str) and cfg_value.upper() == "ALL":
        return ALL_UNIPROT_FEATURE_TYPES
    if isinstance(cfg_value, (list, tuple)):
        return [str(t).upper() for t in cfg_value]
    return ALL_UNIPROT_FEATURE_TYPES

# ----------------------------------------------------------------------------
# Config efetiva a partir de CONFIG
# ----------------------------------------------------------------------------
USER_AGENT = CONFIG["net"]["user_agent"]
DEFAULT_TIMEOUT = (CONFIG["net"]["timeout_connect"], CONFIG["net"]["timeout_read"])

BLAST_DB = CONFIG["blast"]["database"]
BLAST_PROGRAM = "blastp"
BLAST_EXPECT = CONFIG["blast"]["expect"]
BLAST_MAX_ALIGNMENTS = CONFIG["blast"]["hitlist_size"]
BLAST_THROTTLE_SECONDS = CONFIG["blast"]["throttle_seconds"]
MIN_PCT_IDENTITY = CONFIG["blast"]["min_pct_identity"]
MIN_QUERY_COVERAGE = CONFIG["blast"]["min_query_coverage"]
SAME_SPECIES_ONLY = CONFIG["blast"]["same_species_only"]

INTERPRO_PFAM_FOR_PROT = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}?format=json"
PDBe_UNIPROT_MAP = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
UNIPROT_ENTRY_JSON = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_IDMAP_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_IDMAP_STATUS = "https://rest.uniprot.org/idmapping/status/{job}"

EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EUTILS_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EUTILS_ELINK  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
EUTILS_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

NCBI_API_KEY = os.getenv(CONFIG["net"]["ncbi_api_key_env"])
EUTILS_TOOL = "Effatha-RnD"
EUTILS_EMAIL = os.getenv(CONFIG["net"]["eutils_email_env"]) or "lovato.rodrigo@gmail.com"

# ----------------------------------------------------------------------------
# Sessão HTTP com retries
# ----------------------------------------------------------------------------
def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": USER_AGENT, "Accept": "*/*"})
    r = Retry(total=5, backoff_factor=0.5, status_forcelist=(429, 500, 502, 503, 504))
    s.mount("https://", HTTPAdapter(max_retries=r))
    s.mount("http://", HTTPAdapter(max_retries=r))
    return s

SESSION = make_session()

# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------
@dataclass
class PfamDomain:
    accession: str
    start_1based: int
    end_1based: int
    description: str = ""

@dataclass
class Region:
    start_1based: int
    end_1based: int
    source: str             # "UNIPROT_FEATURE", "PFAM", "NCBI_FEATURE"
    type: str               # feature type (ex.: BINDING, ACT_SITE...) ou Pfam ACC
    note: str = ""

@dataclass
class Mapping:
    pdb_id: str
    chain: str
    uniprot_acc: str
    unp_start: int
    unp_end: int
    coverage: float

@dataclass
class Window:
    abs_center_pos_1based: int
    seq: str
    start_in_protein_1based: int
    center_index_in_window_0based: int
    region_source: Optional[str] = None
    region_type: Optional[str] = None
    region_note: Optional[str] = ""

@dataclass
class Variant:
    abs_pos_1based: int
    ref: str
    alt: str
    subject_title: str

# ----------------------------------------------------------------------------
# Globals de mapeamento (para inserir no relatório)
# ----------------------------------------------------------------------------
MAPPING_INFO: Optional[Dict] = None

# ----------------------------------------------------------------------------
# UniProt / InterPro / PDBe utilitários
# ----------------------------------------------------------------------------
def fetch_uniprot_fasta(accession: str) -> str:
    r = SESSION.get(UNIPROT_FASTA_URL.format(acc=accession), timeout=DEFAULT_TIMEOUT)
    r.raise_for_status()
    lines = [ln.strip() for ln in r.text.splitlines() if ln and not ln.startswith(">")]
    return "".join(lines)

def fetch_uniprot_meta(accession: str) -> Dict:
    r = SESSION.get(UNIPROT_ENTRY_JSON.format(acc=accession), timeout=DEFAULT_TIMEOUT)
    r.raise_for_status()
    return r.json()

def fetch_pdb_uniprot_mappings(pdb_id: str) -> List[Mapping]:
    url = PDBe_UNIPROT_MAP.format(pdb_id=pdb_id.lower())
    r = SESSION.get(url, timeout=DEFAULT_TIMEOUT)
    r.raise_for_status()
    data = r.json()
    if not data:
        return []
    root = data.get(pdb_id.lower()) or {}
    uni = root.get("UniProt") or {}
    mappings: List[Mapping] = []
    for acc, info in uni.items():
        for m in info.get("mappings", []):
            chain = m.get("chain_id") or m.get("chain") or "?"
            unp_start = int(m.get("unp_start", 1))
            unp_end = int(m.get("unp_end", unp_start))
            mappings.append(Mapping(pdb_id=pdb_id, chain=chain, uniprot_acc=acc,
                                    unp_start=unp_start, unp_end=unp_end, coverage=0.0))
    return mappings

def choose_best_mapping(maps: List[Mapping], seq_len_by_acc: Dict[str, int], reviewed_by_acc: Dict[str, bool]) -> Mapping:
    best: Optional[Mapping] = None
    for m in maps:
        L = seq_len_by_acc.get(m.uniprot_acc, 1) or 1
        m.coverage = (m.unp_end - m.unp_start + 1) / float(L)
        reviewed = reviewed_by_acc.get(m.uniprot_acc, False)
        if best is None:
            best = m
            best._reviewed = reviewed  # type: ignore
        else:
            if (m.coverage > best.coverage + 1e-9) or (abs(m.coverage - best.coverage) < 1e-9 and reviewed and not getattr(best, "_reviewed", False)):
                best = m
                best._reviewed = reviewed  # type: ignore
    assert best is not None
    return best

def fetch_pfam_domains(accession: str) -> List[PfamDomain]:
    r = SESSION.get(INTERPRO_PFAM_FOR_PROT.format(acc=accession), headers={"Accept": "application/json"}, timeout=DEFAULT_TIMEOUT)
    r.raise_for_status()
    data = r.json()
    out: List[PfamDomain] = []
    for entry in data.get("results", []):
        meta = entry.get("metadata", {})
        pf = meta.get("accession", "")
        desc = meta.get("name", "") or ""
        loc_blocks = []
        if isinstance(entry.get("entry_protein_locations"), list):
            loc_blocks = entry["entry_protein_locations"]
        else:
            for p in entry.get("proteins", []) or []:
                loc_blocks.extend(p.get("entry_protein_locations", []) or [])
        for loc in loc_blocks:
            for frag in loc.get("fragments", []) or []:
                try:
                    s = int(frag.get("start", 0)); e = int(frag.get("end", 0))
                except Exception:
                    continue
                if s > 0 and e >= s:
                    out.append(PfamDomain(pf, s, e, desc))
    return out

# --- UniProt FEATURES (regiões funcionais) ---
def fetch_uniprot_functional_regions(meta: Dict, wanted_types: List[str]) -> List[Region]:
    regions: List[Region] = []
    feats = meta.get("features", []) or []
    for ft in feats:
        ftype = (ft.get("type") or "").upper()
        if ftype not in wanted_types:
            continue
        loc = ft.get("location") or {}
        start = None; end = None
        if "start" in loc and "end" in loc and loc.get("start") and loc.get("end"):
            start = int(loc["start"].get("value", 0)); end = int(loc["end"].get("value", 0))
        elif "position" in loc and loc.get("position"):
            start = end = int(loc["position"].get("value", 0))
        if not start or not end or end < start:
            continue
        notes = []
        desc = ft.get("description")
        if desc:
            notes.append(str(desc))
        lig = ft.get("ligand") or {}
        if isinstance(lig, dict) and (lig.get("name") or lig.get("id") or lig.get("label")):
            name = lig.get("name") or lig.get("label") or "ligand"
            chebi = lig.get("id") or lig.get("databaseId")
            if chebi:
                notes.append(f"ligand={name} ({chebi})")
            else:
                notes.append(f"ligand={name}")
        ligpart = ft.get("ligandPart") or {}
        if isinstance(ligpart, dict) and (ligpart.get("name") or ligpart.get("id")):
            name = ligpart.get("name") or ligpart.get("label") or "part"
            chebi = ligpart.get("id") or ligpart.get("databaseId")
            if chebi:
                notes.append(f"ligandPart={name} ({chebi})")
            else:
                notes.append(f"ligandPart={name}")
        note = "; ".join(notes)
        regions.append(Region(start, end, "UNIPROT_FEATURE", ftype, note))
    return regions

# ----------------------------------------------------------------------------
# NCBI & UniProt mapping helpers
# ----------------------------------------------------------------------------
def _eutils_params(base: Dict[str, str]) -> Dict[str, str]:
    p = dict(base)
    p["tool"] = EUTILS_TOOL
    p["email"] = EUTILS_EMAIL
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p

def _strip_refseq_version(acc: str) -> str:
    # NP_000537.3 -> NP_000537
    if "." in acc:
        return acc.split(".", 1)[0]
    return acc

def _uniprot_idmap_refseq_to_uniprot(acc_list: List[str]) -> Dict[str, str]:
    """
    Usa UniProt ID Mapping (jobs) para mapear RefSeq_Protein -> UniProtKB.
    Retorna dict {input_id: uniprot_acc} (primeiro resultado para cada).
    """
    if not acc_list:
        return {}
    try:
        data = {
            "from": "RefSeq_Protein",
            "to": "UniProtKB",
            "ids": ",".join(acc_list),
        }
        r = SESSION.post(UNIPROT_IDMAP_RUN, data=data, timeout=DEFAULT_TIMEOUT)
        if r.status_code != 200:
            return {}
        job = r.json().get("jobId")
        if not job:
            return {}
        # polling
        for _ in range(30):
            time.sleep(1.0)
            rs = SESSION.get(UNIPROT_IDMAP_STATUS.format(job=job), timeout=DEFAULT_TIMEOUT)
            if rs.status_code != 200:
                continue
            js = rs.json()
            if js.get("jobStatus") == "RUNNING":
                continue
            # resultados inline
            results = (js.get("results") or [])
            out: Dict[str,str] = {}
            for row in results:
                frm = row.get("from")
                to = row.get("to")
                if frm and to and frm not in out:
                    out[frm] = to.split("-", 1)[0]  # garante ACC (sem isoform suffix)
            # também tentar "redirectURL" (algumas versões da API retornam links)
            # mas manteremos simples aqui.
            return out
    except Exception:
        return {}
    return {}

def map_refseq_protein_to_uniprot(refseq_acc: str) -> Optional[str]:
    """
    Estratégia robusta:
      1) UniProt search xref:RefSeq:<acc_versionless>
      2) UniProt search xref:RefSeq:<acc_full>
      3) ID mapping API (com e sem versão)
    """
    candidates: List[str] = []
    for qacc in [_strip_refseq_version(refseq_acc), refseq_acc]:
        try:
            q = f"xref:RefSeq:{qacc}"
            params = {"query": q, "fields": "accession,reviewed", "format": "json", "size": 10}
            r = SESSION.get(UNIPROT_SEARCH_URL, params=params, timeout=DEFAULT_TIMEOUT)
            if r.status_code == 200:
                data = r.json()
                results = data.get("results") or []
                for x in results:
                    acc = x.get("primaryAccession")
                    if acc:
                        candidates.append((acc, (x.get("entryType","").lower()=="reviewed")))
        except Exception:
            pass
        if candidates:
            break

    if not candidates:
        # tenta ID mapping (com e sem versão)
        mp = _uniprot_idmap_refseq_to_uniprot([_strip_refseq_version(refseq_acc), refseq_acc])
        # prioriza mapa do sem versão
        for key in [_strip_refseq_version(refseq_acc), refseq_acc]:
            if key in mp:
                return mp[key]
        # se veio só algum outro, pega o primeiro
        if mp:
            return list(mp.values())[0]
        return None

    # prefere reviewed
    reviewed = [acc for acc,rev in candidates if rev]
    return reviewed[0] if reviewed else candidates[0][0]

# ----------------------------------------------------------------------------
# NCBI: features de proteína (GenPept) e utilitários gene/protein
# ----------------------------------------------------------------------------
def fetch_ncbi_protein_gb_text(acc: str) -> Optional[str]:
    params = _eutils_params({"db": "protein", "id": acc, "rettype": "gp", "retmode": "text"})
    r = SESSION.get(EUTILS_EFETCH, params=params, timeout=DEFAULT_TIMEOUT)
    if r.status_code != 200 or not r.text.strip():
        return None
    return r.text

def fetch_ncbi_protein_features_to_regions(acc: str) -> List[Region]:
    text = fetch_ncbi_protein_gb_text(acc)
    if not text:
        return []
    try:
        rec = SeqIO.read(StringIO(text), "gb")
    except Exception:
        return []
    L = len(rec.seq or "")
    regions: List[Region] = []
    for feat in rec.features or []:
        ftype = (feat.type or "").lower()
        try:
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            if start < 1 or end < start or end > L:
                continue
        except Exception:
            continue

        q = feat.qualifiers or {}
        note_parts: List[str] = []
        if ftype == "site":
            stype = (q.get("site_type",[None])[0] or "").lower()
            if "active" in stype:
                rtype = "ACT_SITE"
            elif "metal" in stype:
                rtype = "BINDING"
                nm = q.get("note",[None])[0]
                if nm:
                    note_parts.append(f"metal={nm}")
            elif "binding" in stype:
                rtype = "BINDING"
            elif "glycosylation" in stype:
                rtype = "CARBOHYD"
            elif "modified" in stype:
                rtype = "MOD_RES"
            else:
                rtype = "SITE"
            if q.get("note"):
                note_parts.append(q["note"][0])
            regions.append(Region(start, end, "NCBI_FEATURE", rtype, "; ".join(note_parts)))
        elif ftype == "bond" or ftype == "disulfide bond":
            rtype = "DISULFID"
            if q.get("note"):
                note_parts.append(q["note"][0])
            regions.append(Region(start, end, "NCBI_FEATURE", rtype, "; ".join(note_parts)))
        elif ftype in ("sig_peptide","signal peptide"):
            regions.append(Region(start, end, "NCBI_FEATURE", "SIGNAL", "; ".join(q.get("note",[]))))
        elif ftype in ("transit_peptide",):
            regions.append(Region(start, end, "NCBI_FEATURE", "TRANSIT", "; ".join(q.get("note",[]))))
        elif ftype in ("mat_peptide",):
            regions.append(Region(start, end, "NCBI_FEATURE", "PEPTIDE", "; ".join(q.get("product",[]))))
        elif ftype in ("propeptide",):
            regions.append(Region(start, end, "NCBI_FEATURE", "PROPEP", "; ".join(q.get("note",[]))))
        elif ftype == "region":
            rname = (q.get("region_name",[None])[0] or q.get("note",[None])[0] or "").lower()
            rtype = "REGION"
            if "zinc finger" in rname:
                rtype = "DOMAIN"; note_parts.append("zinc finger")
            elif "coiled-coil" in rname or "coiled coil" in rname:
                rtype = "COILED"
            elif "transmembrane" in rname:
                rtype = "TRANSMEM"
            elif "dna-binding" in rname or "dna binding" in rname:
                rtype = "DNA_BIND"
            elif "domain" in rname:
                rtype = "DOMAIN"
            if rname:
                note_parts.append(rname)
            regions.append(Region(start, end, "NCBI_FEATURE", rtype, "; ".join(note_parts)))
        elif ftype == "conflict":
            regions.append(Region(start, end, "NCBI_FEATURE", "CONFLICT", "; ".join(q.get("note",[]))))
        elif ftype == "variant":
            regions.append(Region(start, end, "NCBI_FEATURE", "VARIANT", "; ".join(q.get("note",[]))))
        else:
            pass
    return regions

def ncbi_gene_to_refseq_proteins_by_geneid(gene_id: str) -> List[Tuple[str,int]]:
    params = _eutils_params({"dbfrom":"gene","db":"protein","linkname":"gene_protein_refseq","id":gene_id})
    r = SESSION.get(EUTILS_ELINK, params=params, timeout=DEFAULT_TIMEOUT)
    if r.status_code != 200:
        return []
    try:
        root = ET.fromstring(r.text)
        ids = [e.text for e in root.findall(".//LinkSetDb/Link/Id") if e.text]
    except Exception:
        ids = []
    if not ids:
        return []
    params = _eutils_params({"db":"protein","id":",".join(ids)})
    r2 = SESSION.get(EUTILS_ESUMMARY, params=params, timeout=DEFAULT_TIMEOUT)
    if r2.status_code != 200:
        return []
    try:
        root2 = ET.fromstring(r2.text)
        out: List[Tuple[str,int]] = []
        for doc in root2.findall(".//DocSum"):
            acc = None; length = None
            for itm in doc.findall("Item"):
                if itm.get("Name") == "AccessionVersion":
                    acc = itm.text
                if itm.get("Name") == "Length":
                    try:
                        length = int(itm.text or "0")
                    except Exception:
                        length = None
            if acc and length:
                out.append((acc, length))
        return out
    except Exception:
        return []

def ncbi_symbol_to_geneid(symbol: str, taxid: int) -> Optional[str]:
    q = f"{symbol}[Symbol]+AND+{taxid}[TaxID]"
    params = _eutils_params({"db":"gene","term":q,"retmax":"1"})
    r = SESSION.get(EUTILS_ESEARCH, params=params, timeout=DEFAULT_TIMEOUT)
    if r.status_code != 200:
        return None
    try:
        root = ET.fromstring(r.text)
        gid = root.findtext(".//IdList/Id")
        return gid
    except Exception:
        return None

def choose_isoform(refseqs: List[Tuple[str,int]], policy: str = "longest") -> Optional[str]:
    if not refseqs:
        return None
    np = [x for x in refseqs if x[0].startswith("NP_")]
    xp = [x for x in refseqs if x[0].startswith("XP_")]
    group = np if np else xp
    group.sort(key=lambda z: z[1], reverse=True)
    return group[0][0]

# ----------------------------------------------------------------------------
# Regiões → janelas
# ----------------------------------------------------------------------------
def expand_regions(regs: List[Region], L: int, flank_left: int, flank_right: int, merge: bool) -> List[Region]:
    expanded: List[Region] = []
    for r in regs:
        s = max(1, r.start_1based - flank_left)
        e = min(L, r.end_1based + flank_right)
        expanded.append(Region(s, e, r.source, r.type, r.note))
    if not merge or not expanded:
        return expanded
    expanded.sort(key=lambda x: (x.start_1based, x.end_1based))
    merged: List[Region] = [expanded[0]]
    for r in expanded[1:]:
        last = merged[-1]
        if r.start_1based <= last.end_1based + 1 and r.source == last.source and r.type == last.type:
            merged[-1] = Region(last.start_1based, max(last.end_1based, r.end_1based), last.source, last.type, last.note)
        else:
            merged.append(r)
    return merged

def regions_to_windows(sequence: str, regs: List[Region]) -> List[Window]:
    out: List[Window] = []
    for r in regs:
        left = r.start_1based; right = r.end_1based
        subseq = sequence[left - 1: right]
        center0 = (right - left) // 2
        abs_center = left + center0
        out.append(Window(abs_center, subseq, left, center0, r.source, r.type, r.note))
    return out

def make_windows_pfam_mode(sequence: str, domains: List[PfamDomain], flank: int, stride: int) -> List[Window]:
    L = len(sequence)
    win_len = 2 * flank + 1
    windows: List[Window] = []
    for d in domains:
        start, end = d.start_1based, d.end_1based
        pos = start
        while pos <= end:
            left = max(1, pos - flank); right = min(L, pos + flank)
            if right - left + 1 < win_len:
                if left == 1: right = min(L, left + win_len - 1)
                elif right == L: left = max(1, right - win_len + 1)
            subseq = sequence[left - 1: right]
            center0 = pos - left
            windows.append(Window(pos, subseq, left, center0, "PFAM", d.accession, d.description))
            pos += stride
    return windows

# ----------------------------------------------------------------------------
# BLAST
# ----------------------------------------------------------------------------
def qblast_window(
    win: Window,
    *,
    entrez_query: Optional[str],
    min_pct_identity: float,
    min_query_cov: float,
    organism_filter_name: Optional[str] = None,
    hsp_log: bool = False,
) -> List[Variant]:
    kwargs = {
        "program": "blastp",
        "database": BLAST_DB,
        "sequence": win.seq,
        "expect": BLAST_EXPECT,
        "hitlist_size": BLAST_MAX_ALIGNMENTS,
        "format_type": "XML",
    }
    if entrez_query:
        kwargs["entrez_query"] = entrez_query
    if NCBI_API_KEY:
        kwargs["ncbi_api_key"] = NCBI_API_KEY

    try:
        handle = NCBIWWW.qblast(**kwargs)
    except ValueError as e:
        if "Entrez Query" in str(e):
            kwargs.pop("entrez_query", None)
            handle = NCBIWWW.qblast(**kwargs)
        else:
            raise

    vars_out: List[Variant] = []
    for rec in NCBIXML.parse(handle):
        for aln in rec.alignments:
            title = getattr(aln, "title", "")
            if organism_filter_name and organism_filter_name.lower() not in title.lower():
                if hsp_log:
                    print(f"[SKIP-ORG] {title}")
                continue
            for hsp in aln.hsps:
                q_seq = hsp.query or ""
                s_seq = hsp.sbjct or ""
                q_non_gap = sum(1 for c in q_seq if c != "-")
                cov = q_non_gap / max(1, len(win.seq))
                try:
                    pct_id = 100.0 * (hsp.identities / float(hsp.align_length))
                except Exception:
                    pct_id = 0.0
                if cov < min_query_cov or pct_id < min_pct_identity:
                    if hsp_log:
                        print(f"[SKIP] cov={cov:.2f} id={pct_id:.1f}% title={title}")
                    continue
                if hsp_log:
                    print(f"[PASS] cov={cov:.2f} id={pct_id:.1f}% title={title}")
                idx = 0
                for cq, cs in zip(q_seq, s_seq):
                    if cq != "-":
                        abs_pos = win.start_in_protein_1based + idx
                        if cq != cs and cs != "-":
                            vars_out.append(Variant(abs_pos, cq, cs, title))
                        idx += 1
    return vars_out

def qblast_many(
    windows: List[Window],
    *,
    entrez_query: Optional[str],
    min_pct_identity: float,
    min_query_cov: float,
    organism_filter_name: Optional[str] = None,
    hsp_log: bool = False,
    sleep_s: float = BLAST_THROTTLE_SECONDS,
    progress: bool = False,
) -> List[Variant]:
    out: List[Variant] = []
    total = len(windows)
    if progress:
        print(f"[BLAST] Rodando {total} janelas…")
    for i, w in enumerate(windows, 1):
        if progress:
            tag = (w.region_source or "?") + ":" + (w.region_type or "?")
            print(f"  - {i}/{total} pos={w.abs_center_pos_1based} {tag}")
        out.extend(
            qblast_window(
                w,
                entrez_query=entrez_query,
                min_pct_identity=min_pct_identity,
                min_query_cov=min_query_cov,
                organism_filter_name=organism_filter_name,
                hsp_log=hsp_log,
            )
        )
        if i < total:
            time.sleep(sleep_s)
    return out

# ----------------------------------------------------------------------------
# Agregação / Colchetes / CSV
# ----------------------------------------------------------------------------
def variants_by_position(vars: Iterable[Variant]) -> Dict[int, Set[str]]:
    agg: Dict[int, Set[str]] = {}
    for v in vars:
        s = agg.setdefault(v.abs_pos_1based, set())
        s.add(v.ref); s.add(v.alt)
    return agg

def aa_regions_bracketed(sequence: str, regs: List[Region], agg: Dict[int, Set[str]]):
    out = []
    for r in regs:
        chars = []
        for pos in range(r.start_1based, r.end_1based + 1):
            ref = sequence[pos - 1]
            opts = agg.get(pos, set()) or set()
            alts = sorted([a for a in opts if a not in ('-', ref)])
            if alts:
                chars.append('[' + '/'.join([ref] + alts) + ']')
            else:
                chars.append(ref)
        tag = f"{r.source}:{r.type}"
        out.append({"tag": tag, "start": r.start_1based, "end": r.end_1based, "aa": ''.join(chars)})
    return out

def export_variants_csv(vars: List[Variant], path: str) -> None:
    tmp: Dict[int, Dict[str, Set[str]]] = {}
    for v in vars:
        b = tmp.setdefault(v.abs_pos_1based, {"ref": v.ref, "alts": set()})
        if v.alt != v.ref:
            b["alts"].add(v.alt)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["abs_pos", "ref", "alts"])
        for pos in sorted(tmp):
            w.writerow([pos, tmp[pos]["ref"], "/".join(sorted(tmp[pos]["alts"]))])
    print(f"CSV exportado: {path}")

# ----------------------------------------------------------------------------
# DNA/mRNA helpers (NCBI E-utilities) + retrotradução
# ----------------------------------------------------------------------------
def extract_nuccore_ids_from_uniprot_meta(meta: Dict) -> List[str]:
    ids: List[str] = []
    for x in meta.get("uniProtKBCrossReferences", []) or []:
        db = x.get("database")
        if db in ("EMBL", "RefSeq"):
            acc = x.get("id")
            if acc and acc not in ids:
                ids.append(acc)
    return ids

def _parse_fasta_multi(text: str) -> List[Tuple[str, str]]:
    seqs: List[Tuple[str, str]] = []
    header: Optional[str] = None
    buf: List[str] = []
    for ln in text.splitlines():
        if ln.startswith(">"):
            if header is not None:
                seq = "".join(buf).replace(" ", "").replace("\r", "")
                seqs.append((header, seq))
            header = ln[1:].strip()
            buf = []
        else:
            buf.append(ln.strip())
    if header is not None:
        seq = "".join(buf).replace(" ", "").replace("\r", "")
        seqs.append((header, seq))
    return seqs

def _efetch_cds_fasta(nuccore_id: str, rettype: str) -> Optional[str]:
    params = _eutils_params({"db": "nuccore", "id": nuccore_id, "rettype": rettype, "retmode": "text"})
    r = SESSION.get(EUTILS_EFETCH, params=params, timeout=DEFAULT_TIMEOUT)
    if r.status_code != 200:
        return None
    return r.text

def try_fetch_cds_for_uniprot(meta: Dict, protein_seq: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    ids = extract_nuccore_ids_from_uniprot_meta(meta)
    for acc in ids[:5]:
        aa_txt = _efetch_cds_fasta(acc, "fasta_cds_aa")
        if not aa_txt:
            continue
        aa_list = _parse_fasta_multi(aa_txt)
        for hdr, aa_seq in aa_list:
            if aa_seq == protein_seq:
                key = hdr.split(" ")[0]
                na_txt = _efetch_cds_fasta(acc, "fasta_cds_na")
                if not na_txt:
                    continue
                na_list = _parse_fasta_multi(na_txt)
                for hdr2, na_seq in na_list:
                    if hdr2.split(" ")[0] == key and len(na_seq) == 3 * len(aa_seq):
                        return na_seq, acc, key
                for hdr2, na_seq in na_list:
                    if len(na_seq) == 3 * len(aa_seq):
                        return na_seq, acc, hdr2.split(" ")[0]
    return None, None, None

# --- Novo: CDS a partir de RefSeq protein (fallback NCBI-only) ---
def _elink_protein_to_nuccore_ids(refseq_acc: str) -> List[str]:
    # tenta alguns linknames comuns; se falhar, usa o default
    out: List[str] = []
    for linkname in ("protein_nuccore_mrna", "protein_nuccore", None):
        params = _eutils_params({
            "dbfrom": "protein",
            "db": "nuccore",
            "id": refseq_acc,
        })
        if linkname:
            params["linkname"] = linkname
        r = SESSION.get(EUTILS_ELINK, params=params, timeout=DEFAULT_TIMEOUT)
        if r.status_code != 200:
            continue
        try:
            root = ET.fromstring(r.text)
            ids = [e.text for e in root.findall(".//LinkSetDb/Link/Id") if e.text]
            out.extend(ids)
        except Exception:
            continue
    # remover duplicatas preservando ordem
    seen=set(); uniq=[]
    for x in out:
        if x not in seen:
            seen.add(x); uniq.append(x)
    return uniq

def try_fetch_cds_for_refseq_protein(refseq_acc: str, protein_seq: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    nuccores = _elink_protein_to_nuccore_ids(refseq_acc)
    for acc in nuccores[:10]:
        aa_txt = _efetch_cds_fasta(acc, "fasta_cds_aa")
        if not aa_txt:
            continue
        aa_list = _parse_fasta_multi(aa_txt)
        hit_keys = []
        for hdr, aa_seq in aa_list:
            if aa_seq == protein_seq:
                hit_keys.append(hdr.split(" ")[0])
        if not hit_keys:
            continue
        na_txt = _efetch_cds_fasta(acc, "fasta_cds_na")
        if not na_txt:
            continue
        na_list = _parse_fasta_multi(na_txt)
        for hdr2, na_seq in na_list:
            key2 = hdr2.split(" ")[0]
            if key2 in hit_keys and len(na_seq) == 3 * len(protein_seq):
                return na_seq, acc, key2
    return None, None, None

_CODON_PREF = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT",
    "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
    "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT",
    "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTG",
    "*": "TAA",
}

def naive_reverse_translate(protein_seq: str) -> str:
    return "".join(_CODON_PREF.get(aa, "NNN") for aa in protein_seq)

def nt_segments_for_regions(regs: List[Region], cds_na_seq: str) -> List[Tuple[str, int, int, str, str, str]]:
    """
    Retorna lista de tuplas (tag, start, end, dna, mrna, note),
    onde tag = "SOURCE:TYPE" e note carrega a descrição/anotação da região.
    """
    segs: List[Tuple[str, int, int, str, str, str]] = []
    for r in regs:
        s = max(0, (r.start_1based - 1) * 3)
        e = min(len(cds_na_seq), r.end_1based * 3)
        if s >= e:
            continue
        dna = cds_na_seq[s:e]
        mrna = dna.replace("T", "U").replace("t", "u")
        note = r.note or ""
        segs.append((f"{r.source}:{r.type}", r.start_1based, r.end_1based, dna, mrna, note))
    return segs

# ----------------------------------------------------------------------------
# API watch (placeholder)
# ----------------------------------------------------------------------------
def check_api_updates(config: Dict, format_vars: Optional[Dict] = None) -> Optional[str]:
    return None

# ----------------------------------------------------------------------------
# Relatório .txt
# ----------------------------------------------------------------------------
def write_report_txt(path: str, *, context: Dict) -> None:
    """Gera um relatório .txt consolidando config, resumo, regiões e sequências, com notas das regiões no cabeçalho de cada bloco."""
    lines: List[str] = []
    lines.append("# Relatório de Execução — Pipeline funcional")
    lines.append(f"Data/hora: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("## Descrição do processo")
    lines.append("1) Identificação da proteína (PDB/UniProt/NCBI) e espécie.")
    lines.append("2) Coleta de regiões funcionais (UniProt / Pfam / NCBI).")
    lines.append("3) Geração de janelas pelas regiões ± flancos.")
    lines.append("4) BLASTP intra-espécie por janela (opcional).")
    lines.append("5) Agregação por posição absoluta e anotação em colchetes (ref primeiro).")
    lines.append("6) Mapeamento para DNA/mRNA real (CDS) quando possível; fallback retrotradução.")
    lines.append("")

    # Config (resumo seguro)
    lines.append("## Config usada (resumo)")
    try:
        lines.append(json.dumps({
            "input.source": CONFIG["input"]["source"],
            "regions.mode": CONFIG["regions"]["mode"],
            "regions.include_ncbi_protein_features": CONFIG["regions"].get("include_ncbi_protein_features", False),
            "blast.enable": CONFIG["blast"]["enable"],
            "variation_layers": {k:v["enable"] for k,v in CONFIG.get("external_layers",{}).items()},
            "filters.enabled": CONFIG["variation_filters"]["enabled"],
        }, indent=2))
    except Exception:
        lines.append(json.dumps(CONFIG, indent=2))
    lines.append("")

    # Fonte & Mapeamento (quando existir)
    mapping = context.get("mapping") or {}
    if mapping:
        lines.append("## Fonte & Mapeamento")
        for k in ["input_source","input_id","resolved_gene","resolved_refseq","mapped_uniprot","notes"]:
            if mapping.get(k) is not None:
                lines.append(f"- {k}: {mapping.get(k)}")
        lines.append("")

    # Proteína
    prot = context.get("protein", {})
    lines.append("## Proteína")
    lines.append(f"Acesso: {prot.get('accession','?')}")
    lines.append(f"Nome: {prot.get('name','')}")
    lines.append(f"Espécie: {prot.get('organism','?')} (taxid={prot.get('taxid','?')})")
    lines.append(f"Tamanho (aa): {prot.get('length','?')}")
    lines.append("")

    # Regiões
    lines.append("## Regiões-alvo (após expansão)")
    for r in context.get("regions", []):
        note = f" — {r['note']}" if r.get('note') else ""
        lines.append(f"- {r['source']}:{r['type']} {r['start']}-{r['end']}{note}")
    lines.append("")

    # BLAST resumo
    bl = context.get("blast", {})
    lines.append("## BLAST")
    lines.append(f"Habilitado: {'SIM' if CONFIG['blast']['enable'] else 'NÃO'}; Janelas processadas: {bl.get('windows',0)}; Banco: {BLAST_DB}; Filtro: {bl.get('filter','None')}")
    lines.append(f"Identidade mínima: {MIN_PCT_IDENTITY}%; Cobertura mínima: {MIN_QUERY_COVERAGE}")
    lines.append(f"Variações (posições únicas) — BLAST: {bl.get('variant_positions_blast', 0)}; UniProt VARIANT: {bl.get('variant_positions_uniprot', 0)}; Proteins Variation: {bl.get('variant_positions_proteins', 0)}; Filtradas: {bl.get('variant_positions_filtered', 0)}")
    lines.append("")

    # DNA/mRNA + AA por região (com nota/descrição junto ao cabeçalho)
    lines.append("## DNA/mRNA por região (completo) + AA anotado")
    aa_index = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions", [])}
    if not aa_index and context.get("aa_regions_obs"):
        aa_index = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions_obs", [])}

    for seg in context.get("nt_segments", []):
        key = (seg["tag"], seg["start"], seg["end"])
        aa_line = aa_index.get(key, "")
        header = f"{seg['tag']} {seg['start']}-{seg['end']}"
        if seg.get("note"):
            header += f" — {seg['note']}"
        lines.append(header)
        if aa_line:
            lines.append(f"AA  : {aa_line}")
        lines.append(f"DNA : {seg['dna']}")
        lines.append(f"mRNA: {seg['mrna']}")
        lines.append("")

    # Sequência AA completa
    if prot.get("sequence"):
        lines.append("## Sequência de aminoácidos (completa)")
        lines.append(prot["sequence"])
        lines.append("")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Relatório salvo em: {path}")

# ----------------------------------------------------------------------------
# Orquestração — pipeline com UniProt
# ----------------------------------------------------------------------------
def build_regions(sequence: str, accession: str, meta: Dict) -> List[Region]:
    mode = CONFIG["regions"]["mode"]
    regs: List[Region] = []

    if mode in ("uniprot_features", "both"):
        ft_types = _resolve_feature_types(CONFIG["regions"]["feature_types"])
        regs.extend(fetch_uniprot_functional_regions(meta, ft_types))

    if mode in ("pfam", "both"):
        pf = fetch_pfam_domains(accession)
        for d in pf:
            regs.append(Region(d.start_1based, d.end_1based, "PFAM", d.accession, d.description))

    if CONFIG["regions"].get("include_ncbi_protein_features", False):
        refseqs = extract_refseq_proteins_from_uniprot_meta(meta, prefer_np=True)
        added = 0
        for rs in refseqs[:3]:
            ncbi_regs = fetch_ncbi_protein_features_to_regions(rs)
            regs.extend(ncbi_regs)
            added += len(ncbi_regs)
        if added:
            print(f"[INFO] NCBI features agregadas: {added} regiões (via RefSeq cross-ref).")

    L = len(sequence)
    regs = expand_regions(
        regs,
        L,
        CONFIG["regions"]["flank_left"],
        CONFIG["regions"]["flank_right"],
        CONFIG["regions"]["merge_overlaps"],
    )
    return regs

def extract_refseq_proteins_from_uniprot_meta(meta: Dict, prefer_np: bool = True) -> List[str]:
    ids: List[str] = []
    for x in meta.get("uniProtKBCrossReferences", []) or []:
        if x.get("database") == "RefSeq":
            acc = x.get("id")
            if acc and (acc.startswith("NP_") or acc.startswith("XP_")):
                ids.append(acc)
    if prefer_np:
        nps = [a for a in ids if a.startswith("NP_")]
        xps = [a for a in ids if a.startswith("XP_")]
        return nps + xps
    return ids

def _core_pipeline_from_uniprot(accession: str, *, mapping_info: Optional[Dict]=None) -> None:
    # Sequência + metadados UniProt
    seq = fetch_uniprot_fasta(accession)
    print(f"Tamanho da proteína: {len(seq)} aa")
    meta = fetch_uniprot_meta(accession)

    # Espécie / nome
    org = meta.get("organism", {})
    sci_name = org.get("scientificName") or ""
    taxid = None
    try:
        taxid = int(org.get("taxonId"))
    except Exception:
        pass
    prot_name = (meta.get("proteinDescription", {}).get("recommendedName", {}) or {}).get("fullName", {}).get("value", "")
    print(f"Espécie: {sci_name} (taxid={taxid})")

    # Regiões (UniProt/Pfam/NCBI, conforme sua build_regions)
    regs = build_regions(seq, accession, meta)
    print(f"Regiões-alvo (após expansão): {len(regs)}")
    for r in regs[:10]:
        extra = f" — {r.note}" if r.note else ""
        print(f"  - {r.source}:{r.type} {r.start_1based}-{r.end_1based}{extra}")
    if len(regs) > 10:
        print(f"  … (+{len(regs)-10} regiões)")

    windows = regions_to_windows(seq, regs) if regs else []
    print(f"Janelas: {len(windows)}")

    # BLAST opcional
    vars_: List[Variant] = []
    if CONFIG["blast"].get("enable", False):
        entrez = f"txid{taxid}[ORGN]" if (SAME_SPECIES_ONLY and taxid) else None
        print(f"BLAST: {len(windows)} janelas; filtro={entrez}")
        vars_ = qblast_many(
            windows,
            entrez_query=entrez,
            min_pct_identity=MIN_PCT_IDENTITY,
            min_query_cov=MIN_QUERY_COVERAGE,
            organism_filter_name=sci_name if SAME_SPECIES_ONLY else None,
            hsp_log=False,
            progress=CONFIG["output"]["blast_progress"],
        )
    else:
        print("BLAST desabilitado por configuração.")

    agg = variants_by_position(vars_)

    # AA anotado + NT por região (com nota)
    regs_for_nt = regs if regs else [Region(1, len(seq), "FULL", "FULL")]
    aa_by_region = aa_regions_bracketed(seq, regs_for_nt, agg)

    nt_segments: List[Tuple[str, int, int, str, str, str]] = []
    if CONFIG["cds_mapping"]["enable"]:
        cds_na, nuccore_src, header_key = try_fetch_cds_for_uniprot(meta, seq)
        if cds_na and len(cds_na) == 3 * len(seq):
            print(f"CDS encontrado em {nuccore_src} (nt={len(cds_na)}).")
            nt_segments = nt_segments_for_regions(regs_for_nt, cds_na)
        else:
            print("Não foi possível mapear CDS real; usando retrotradução (código genético padrão).")
            dna_full = naive_reverse_translate(seq)
            nt_segments = nt_segments_for_regions(regs_for_nt, dna_full)
    else:
        dna_full = naive_reverse_translate(seq)
        nt_segments = nt_segments_for_regions(regs_for_nt, dna_full)

    # Impressão final — inclui nota no cabeçalho
    print("\n## DNA/mRNA por região (completo) + AA anotado")
    aa_index = {(a["tag"], a["start"], a["end"]): a["aa"] for a in aa_by_region}
    for tag, ds, de, dna_seg, mrna_seg, note in nt_segments:
        header = f"{tag} {ds}-{de}"
        if note:
            header += f" — {note}"
        print(header)
        aa_line = aa_index.get((tag, ds, de), "")
        if aa_line:
            print(f"AA  : {aa_line}")
        print(f"DNA : {dna_seg}")
        print(f"mRNA: {mrna_seg}")

    # Contexto para UI/artefatos — inclui note nos segmentos
    context = {
        "mapping": mapping_info or {},
        "protein": {
            "accession": accession,
            "name": prot_name,
            "organism": sci_name,
            "taxid": taxid,
            "length": len(seq),
            "sequence": seq,
        },
        "regions": [
            {"source": r.source, "type": r.type, "start": r.start_1based, "end": r.end_1based, "note": r.note}
            for r in regs
        ],
        "blast": {
            "windows": len(windows),
            "filter": (f"txid{taxid}[ORGN]" if (CONFIG['blast'].get('enable', False) and taxid) else "None"),
            "variant_positions": len({p for p in agg.keys()}),

            # chaves usadas pela UI (preencha se você agrega externas em outro lugar)
            "variant_positions_blast": len({p for p in agg.keys()}),
            "variant_positions_uniprot": 0,
            "variant_positions_proteins": 0,
            "variant_positions_filtered": 0,
        },
        "layers": {
            "uniprot_variant_enabled": bool(CONFIG.get("external_layers",{}).get("uniprot_variant",{}).get("enable", False)),
            "proteins_variation_enabled": bool(CONFIG.get("external_layers",{}).get("variation_api",{}).get("enable", False)),
            "blast_enabled": bool(CONFIG.get("blast",{}).get("enable", False)),
        },
        "aa_regions": aa_by_region,
        "aa_regions_obs": [],
        "aa_regions_filtered": [],
        "nt_segments": [
            {"tag": tag, "start": ds, "end": de, "dna": dna, "mrna": mrna, "note": note}
            for (tag, ds, de, dna, mrna, note) in nt_segments
        ],
    }

    report_path = CONFIG["output"]["report_txt"]
    if report_path:
        write_report_txt(report_path, context=context)

    arts = CONFIG["output"].get("artifacts_dir")
    if arts:
        try:
            with open(os.path.join(arts, "context_summary.json"), "w", encoding="utf-8") as fh:
                json.dump(context, fh, ensure_ascii=False, indent=2)
        except Exception:
            pass

# ----------------------------------------------------------------------------
# Fallback NCBI-only (sem UniProt)
# ----------------------------------------------------------------------------
def _core_pipeline_from_refseq_only(refseq_acc: str, *, mapping_info: Optional[Dict]=None) -> None:
    print("[WARN] Executando em modo NCBI-only (sem mapeamento UniProt).")

    # GenPept
    gp_txt = fetch_ncbi_protein_gb_text(refseq_acc)
    if not gp_txt:
        print("[ERRO] Não foi possível obter GenPept para a proteína RefSeq.")
        return
    rec = SeqIO.read(StringIO(gp_txt), "gb")
    seq = str(rec.seq)
    print(f"Tamanho da proteína: {len(seq)} aa")
    org_name = rec.annotations.get("organism") or "?"

    # TaxID via ESummary
    taxid = None
    try:
        params = _eutils_params({"db":"protein","id":refseq_acc})
        r = SESSION.get(EUTILS_ESUMMARY, params=params, timeout=DEFAULT_TIMEOUT)
        if r.status_code == 200:
            root = ET.fromstring(r.text)
            for doc in root.findall(".//DocSum"):
                for itm in doc.findall("Item"):
                    if itm.get("Name") == "TaxId":
                        taxid = int(itm.text or "0")
                        break
    except Exception:
        pass
    print(f"Espécie: {org_name} (taxid={taxid})")

    # Regiões a partir das features NCBI
    ncbi_regs = fetch_ncbi_protein_features_to_regions(refseq_acc)
    regs = expand_regions(
        ncbi_regs, len(seq),
        CONFIG["regions"]["flank_left"],
        CONFIG["regions"]["flank_right"],
        CONFIG["regions"]["merge_overlaps"],
    )
    print(f"Regiões-alvo (NCBI features, após expansão): {len(regs)}")
    for r in regs[:10]:
        extra = f" — {r.note}" if r.note else ""
        print(f"  - {r.source}:{r.type} {r.start_1based}-{r.end_1based}{extra}")
    if len(regs) > 10:
        print(f"  … (+{len(regs)-10} regiões)")

    windows = regions_to_windows(seq, regs) if regs else []
    print(f"Janelas: {len(windows)}")

    # BLAST opcional
    vars_: List[Variant] = []
    if CONFIG["blast"].get("enable", False):
        entrez = f"txid{taxid}[ORGN]" if (SAME_SPECIES_ONLY and taxid) else None
        print(f"BLAST: {len(windows)} janelas; filtro={entrez}")
        vars_ = qblast_many(
            windows,
            entrez_query=entrez,
            min_pct_identity=MIN_PCT_IDENTITY,
            min_query_cov=MIN_QUERY_COVERAGE,
            organism_filter_name=org_name if SAME_SPECIES_ONLY else None,
            hsp_log=False,
            progress=CONFIG["output"]["blast_progress"],
        )
    else:
        print("BLAST desabilitado por configuração.")

    agg = variants_by_position(vars_)

    # AA anotado + NT por região (com nota)
    regs_for_nt = regs if regs else [Region(1, len(seq), "FULL", "FULL")]
    aa_by_region = aa_regions_bracketed(seq, regs_for_nt, agg)

    nt_segments: List[Tuple[str, int, int, str, str, str]] = []
    if CONFIG["cds_mapping"]["enable"]:
        cds_na, nuccore_src, header_key = try_fetch_cds_for_refseq_protein(refseq_acc, seq)
        if cds_na and len(cds_na) == 3 * len(seq):
            print(f"CDS encontrado em {nuccore_src} (nt={len(cds_na)}).")
            nt_segments = nt_segments_for_regions(regs_for_nt, cds_na)
        else:
            print("Não foi possível mapear CDS real; usando retrotradução (código genético padrão).")
            dna_full = naive_reverse_translate(seq)
            nt_segments = nt_segments_for_regions(regs_for_nt, dna_full)
    else:
        dna_full = naive_reverse_translate(seq)
        nt_segments = nt_segments_for_regions(regs_for_nt, dna_full)

    # Impressão final — inclui nota no cabeçalho
    print("\n## DNA/mRNA por região (completo) + AA anotado")
    aa_index = {(a["tag"], a["start"], a["end"]): a["aa"] for a in aa_by_region}
    for tag, ds, de, dna_seg, mrna_seg, note in nt_segments:
        header = f"{tag} {ds}-{de}"
        if note:
            header += f" — {note}"
        print(header)
        aa_line = aa_index.get((tag, ds, de), "")
        if aa_line:
            print(f"AA  : {aa_line}")
        print(f"DNA : {dna_seg}")
        print(f"mRNA: {mrna_seg}")

    # Contexto para UI/artefatos
    context = {
        "mapping": mapping_info or {},
        "protein": {
            "accession": refseq_acc,   # no modo NCBI-only mostramos o RefSeq
            "name": rec.name or "",
            "organism": org_name,
            "taxid": taxid,
            "length": len(seq),
            "sequence": seq,
        },
        "regions": [
            {"source": r.source, "type": r.type, "start": r.start_1based, "end": r.end_1based, "note": r.note}
            for r in regs
        ],
        "blast": {
            "windows": len(windows),
            "filter": (f"txid{taxid}[ORGN]" if (CONFIG['blast'].get('enable', False) and taxid) else "None"),
            "variant_positions": len({p for p in agg.keys()}),

            "variant_positions_blast": len({p for p in agg.keys()}),
            "variant_positions_uniprot": 0,
            "variant_positions_proteins": 0,
            "variant_positions_filtered": 0,
        },
        "layers": {
            "uniprot_variant_enabled": bool(CONFIG.get("external_layers",{}).get("uniprot_variant",{}).get("enable", False)),
            "proteins_variation_enabled": bool(CONFIG.get("external_layers",{}).get("variation_api",{}).get("enable", False)),
            "blast_enabled": bool(CONFIG.get("blast",{}).get("enable", False)),
        },
        "aa_regions": aa_by_region,
        "aa_regions_obs": [],
        "aa_regions_filtered": [],
        "nt_segments": [
            {"tag": tag, "start": ds, "end": de, "dna": dna, "mrna": mrna, "note": note}
            for (tag, ds, de, dna, mrna, note) in nt_segments
        ],
    }

    report_path = CONFIG["output"]["report_txt"]
    if report_path:
        write_report_txt(report_path, context=context)

    arts = CONFIG["output"].get("artifacts_dir")
    if arts:
        try:
            with open(os.path.join(arts, "context_summary.json"), "w", encoding="utf-8") as fh:
                json.dump(context, fh, ensure_ascii=False, indent=2)
        except Exception:
            pass

# ----------------------------------------------------------------------------
# Entradas
# ----------------------------------------------------------------------------
def run_from_uniprot(accession: str) -> None:
    print(f"[INFO] Usando UniProt {accession}…")
    _core_pipeline_from_uniprot(accession, mapping_info=MAPPING_INFO)

def run_from_pdb(pdb_id: str) -> None:
    pdb_id = pdb_id.strip().lower()
    print(f"Descobrindo UniProt/espécie para PDB {pdb_id} via SIFTS…")
    maps = fetch_pdb_uniprot_mappings(pdb_id)
    if not maps:
        print("Nenhum mapeamento UniProt encontrado para este PDB.")
        return
    seq_cache: Dict[str, str] = {}
    rev_cache: Dict[str, bool] = {}
    for m in maps:
        if m.uniprot_acc not in seq_cache:
            seq_cache[m.uniprot_acc] = fetch_uniprot_fasta(m.uniprot_acc)
        if m.uniprot_acc not in rev_cache:
            meta_tmp = fetch_uniprot_meta(m.uniprot_acc)
            rev_cache[m.uniprot_acc] = (meta_tmp.get("entryType", "").lower().find("reviewed") >= 0)
    seq_len_by_acc = {acc: len(seq_cache[acc]) for acc in seq_cache}
    best = choose_best_mapping(maps, seq_len_by_acc, rev_cache)
    acc = best.uniprot_acc
    print(f"Selecionado: UniProt {acc} (cadeia {best.chain})")
    _core_pipeline_from_uniprot(acc, mapping_info=MAPPING_INFO)

def run_from_ncbi_protein(refseq_acc: str) -> None:
    global MAPPING_INFO  # <-- declare no topo da função
    refseq_acc = refseq_acc.strip()
    print(f"[INFO] Usando NCBI Protein {refseq_acc}…")

    unp = map_refseq_protein_to_uniprot(refseq_acc)
    if not unp:
        print("[WARN] Não foi possível mapear RefSeq→UniProt. Seguindo em modo NCBI-only.")
        MAPPING_INFO = {
            "input_source": "ncbi_protein",
            "input_id": refseq_acc,
            "resolved_refseq": refseq_acc,
            "mapped_uniprot": None,
            "notes": "NCBI-only fallback (sem UniProt)"
        }
        _core_pipeline_from_refseq_only(refseq_acc, mapping_info=MAPPING_INFO)
        return

    print(f"[INFO] Mapeado para UniProt {unp}")
    MAPPING_INFO = {
        "input_source": "ncbi_protein",
        "input_id": refseq_acc,
        "resolved_refseq": refseq_acc,
        "mapped_uniprot": unp,
        "notes": "mapeado via UniProt (xref/idmapping)"
    }
    _core_pipeline_from_uniprot(unp, mapping_info=MAPPING_INFO)

def run_from_ncbi_gene(gene_cfg: Dict) -> None:
    global MAPPING_INFO  # <-- declare no topo da função
    print(f"[INFO] Usando NCBI Gene…")

    gid = None
    if (gene_cfg.get("id_type") == "entrez") and gene_cfg.get("id"):
        gid = gene_cfg.get("id")
    elif (gene_cfg.get("id_type") == "symbol") and gene_cfg.get("symbol") and gene_cfg.get("taxid"):
        gid = ncbi_symbol_to_geneid(gene_cfg["symbol"], int(gene_cfg["taxid"]))
    if not gid:
        print("[ERRO] Gene não resolvido (informe GeneID ou Símbolo+TaxID).")
        return

    refseqs = ncbi_gene_to_refseq_proteins_by_geneid(gid)
    if not refseqs:
        print("[ERRO] Nenhuma proteína RefSeq associada ao gene.")
        return

    acc = choose_isoform(refseqs, policy=gene_cfg.get("isoform_policy","longest"))
    if not acc:
        print("[ERRO] Falha ao escolher isoforma.")
        return
    print(f"[INFO] Isoforma escolhida: {acc}")

    unp = map_refseq_protein_to_uniprot(acc)
    if not unp:
        print("[WARN] Não foi possível mapear isoforma para UniProt. Seguindo em modo NCBI-only.")
        MAPPING_INFO = {
            "input_source": "ncbi_gene",
            "input_id": (gene_cfg.get("id") or gene_cfg.get("symbol")),
            "resolved_gene": f"GeneID:{gid}",
            "resolved_refseq": acc,
            "mapped_uniprot": None,
            "notes": "gene -> refseq protein (policy) -> NCBI-only fallback"
        }
        _core_pipeline_from_refseq_only(acc, mapping_info=MAPPING_INFO)
        return

    print(f"[INFO] Mapeado para UniProt {unp}")
    MAPPING_INFO = {
        "input_source": "ncbi_gene",
            "input_id": (gene_cfg.get("id") or gene_cfg.get("symbol")),
            "resolved_gene": f"GeneID:{gid}",
            "resolved_refseq": acc,
            "mapped_uniprot": unp,
            "notes": "gene -> refseq protein (policy) -> uniprot (xref/idmapping)"
    }
    _core_pipeline_from_uniprot(unp, mapping_info=MAPPING_INFO)

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    p = argparse.ArgumentParser()
    p.add_argument("--pdb", help="PDB ID (ex.: 5xnl)")
    p.add_argument("--uniprot", help="UniProt ACC (ex.: P02957)")
    p.add_argument("--ncbi_protein", help="NCBI RefSeq protein (ex.: NP_000537.3)")
    p.add_argument("--gene_id", help="NCBI GeneID (ex.: 7157)")
    p.add_argument("--gene_symbol", help="Gene symbol (ex.: TP53)")
    p.add_argument("--gene_taxid", help="Gene TaxID (ex.: 9606)")
    args = p.parse_args()

    source = CONFIG["input"]["source"].lower()
    if args.pdb:
        source = "pdb"; CONFIG["input"]["pdb_id"] = args.pdb
    if args.uniprot:
        source = "uniprot"; CONFIG["input"]["uniprot_acc"] = args.uniprot
    if args.ncbi_protein:
        source = "ncbi_protein"; CONFIG["input"]["ncbi_protein_acc"] = args.ncbi_protein
    if args.gene_id or args.gene_symbol:
        source = "ncbi_gene"
        if args.gene_id:
            CONFIG["input"]["gene"]["id_type"] = "entrez"
            CONFIG["input"]["gene"]["id"] = args.gene_id
        else:
            CONFIG["input"]["gene"]["id_type"] = "symbol"
            CONFIG["input"]["gene"]["symbol"] = args.gene_symbol or ""
            CONFIG["input"]["gene"]["taxid"] = int(args.gene_taxid or "9606")

    if source == "pdb":
        run_from_pdb(CONFIG["input"]["pdb_id"])
    elif source == "uniprot":
        run_from_uniprot(CONFIG["input"]["uniprot_acc"])
    elif source == "ncbi_protein":
        run_from_ncbi_protein(CONFIG["input"]["ncbi_protein_acc"])
    elif source == "ncbi_gene":
        run_from_ncbi_gene(CONFIG["input"]["gene"])
    else:
        raise ValueError("CONFIG['input']['source'] deve ser 'pdb' | 'uniprot' | 'ncbi_protein' | 'ncbi_gene'.")

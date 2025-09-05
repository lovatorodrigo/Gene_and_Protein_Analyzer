#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Effatha — Gene & Protein Analyzer (com regiões funcionais em todas as entradas)
- Normaliza qualquer entrada (UniProt, PDB, NCBI Protein, NCBI Gene) para um alvo UniProt
- Extrai features reais do UniProt como regiões/sítios (com flancos configuráveis)
- Mapeia CDS real (GenBank/nuccore) sem retrotradução (por /protein_id ou /translation)
- BLASTp/BLASTn por região → sintaxe Effatha [A/L/E] (AA) e [A/C/G/T] (NT)
- Saídas: runs/report.txt, runs/regions.csv, runs/context_summary.json (+ opcional variants_blast.csv)

Requisitos:
  pip install biopython requests
  export NCBI_EMAIL="seu_email@dominio"
  export NCBI_API_KEY="sua_chave"   # opcional, mas recomendado
"""

from __future__ import annotations

import os, io, re, csv, json, time
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any

import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

# =========================
# Config
# =========================
CONFIG: Dict[str, Any] = {
    "input": {
        "source": "ncbi_protein",  # "uniprot" | "pdb" | "ncbi_protein" | "ncbi_gene"
        "uniprot_acc": "P00533",
        "pdb_id": "5XNL",
        "ncbi_protein_acc": "AML61188.1",
        "gene": {
            "id_type": "entrez",  # "entrez" | "symbol"
            "id": "1956",             # ex.: "7157"
            "symbol": "",         # ex.: "TP53"
            "taxid": 9606,
            "isoform_policy": "longest",  # usado apenas para entrada "ncbi_gene"
        },
    },

    # Regiões/sítios funcionais (a partir de features reais do UniProt)
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
        "point_flank": 5,         # flancos SOMENTE em features pontuais (±5 aa)
        "default_min_len": 6,     # descarta regiões curtas após aplicar flancos
        "merge_overlaps": True,   # mesclar sobreposições do MESMO tipo
        "add_full": True,         # adiciona FULL (1..len)
        # NOVO: se o UniProt vier vazio (ou só FULL), tenta features do GenPept automaticamente
        "fallback_genpept_if_uniprot_featureless": True
    },

    "blast": {
        "enable": True,  # liga/desliga BLAST (AA + NT)

        # ===== Controles de espécie/preview =====
        "same_species_only": True,         # aplica filtro de mesma espécie em AA e NT (via taxid)
        "log_offspecies_preview": False,   # roda um preview sem filtro (somente LOG; não entra nos colchetes)

        # ===== Preset para MÁXIMO de colchetes no MESMO organismo =====
        "protein": {
            # roda em múltiplos bancos e une as variações (mesmo organismo via entrez_query)
            "dbs": ["refseq_protein", "nr"],
            # "db" single é mantido como fallback (não usado quando "dbs" presente)
            "db": "refseq_protein",

            "hitlist_size": 200,
            "expect": 1e-5,

            # thresholds permissivos para capturar isoformas/variantes intra-espécie
            "min_identity": 0.85,          # proporção
            "min_query_coverage": 0.90,

            # será setado dinamicamente quando tivermos taxid (ex.: "txid9606[ORGN]")
            "entrez_query": None,
        },
        "nt": {
            "dna_db": "nt",           # DNA
            "rna_db": "refseq_rna",   # mRNA
            "hitlist_size": 25,
            "expect": 1e-10,
            "min_identity": 0.98,
            "min_query_coverage": 0.95,
            "megablast": True,
            # NOVO: filtro por espécie também no NT
            "entrez_query": None
        }
    },

    "cds_mapping": {"enable": True},

    "output": {
        "artifacts_dir": "runs",
        "report_txt": "runs/report.txt",
        "export_regions_csv": "runs/regions.csv",
        "export_csv": None,  # "runs/variants_blast.csv"
        "blast_progress": True,

        # ===== Auditoria de BLAST =====
        "log_hits": True,                           # imprime hits aceitos no console
        "blast_hits_csv": "runs/blast_hits.csv",    # salva auditoria de hits em CSV
    },

    # (opcional) alternativa para setar credenciais aqui em vez de ENV:
    "ncbi": {
        # "email": "seu_email@dominio",
        # "api_key": "sua_chave_aqui"
    }
}

# =========================
# NCBI / HTTP
# =========================
NCBI_EMAIL   = os.getenv("NCBI_EMAIL", "").strip()
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()
NCBI_EMAIL   = CONFIG.get("ncbi", {}).get("email", NCBI_EMAIL).strip()
NCBI_API_KEY = CONFIG.get("ncbi", {}).get("api_key", NCBI_API_KEY).strip()

if NCBI_EMAIL:
    Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Effatha-GPA/3.0 (+https://effatha)"})


def _http_retry(method: str, url: str, **kwargs):
    """GET/POST com backoff (429/5xx)."""
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
            if i == max_tries - 1:
                raise
            time.sleep((backoff ** i) + (0.1 * i))
        except requests.RequestException:
            if i == max_tries - 1:
                raise
            time.sleep((backoff ** i) + (0.1 * i))
    raise RuntimeError("HTTP retry failed")

# =========================
# Estruturas
# =========================
@dataclass
class Region:
    start_1based: int
    end_1based: int
    tag: str
    note: str = ""

    @property
    def length(self) -> int:
        return self.end_1based - self.start_1based + 1


@dataclass
class PDBUniProtMap:
    uniprot_acc: str
    chain: str
    coverage: float     # 0..1

# =========================
# Auditoria de BLAST (logs/CSV) + alvo esperado
# =========================
_BLASTP_AUDIT: List[Dict[str, Any]] = []
_BLASTN_AUDIT: List[Dict[str, Any]] = []

_EXPECTED_TAXID: Optional[int] = None
_EXPECTED_SPECIES: Optional[str] = None

def _species_from_def(s: str) -> str:
    m = re.search(r"\[([^\[\]]+)\]\s*$", s or "")
    return m.group(1) if m else ""

def _audit_hit(kind: str, db: str, region_tag: str, accession: str,
               hit_def: str, pid: float, cov: float, hsp_query_len: int):
    """
    kind: "protein" ou "nucleotide"
    pid, cov: 0..1
    - extrai espécie do defline (quando presente no padrão "... [Species]")
    - marca same_species se coincidir com o alvo atual (quando conhecido)
    """
    hit_species = _species_from_def(hit_def)
    same_species = None
    if _EXPECTED_SPECIES:
        same_species = (hit_species == _EXPECTED_SPECIES)

    rec = {
        "kind": kind,
        "db": db,
        "region": region_tag,
        "accession": accession,
        "percent_identity": round(pid * 100.0, 1),
        "query_coverage": round(cov * 100.0, 1),
        "hsp_query_len": int(hsp_query_len),
        "definition": (hit_def or "")[:200],
        "species": hit_species,
        "same_species": same_species
    }
    if kind == "protein":
        _BLASTP_AUDIT.append(rec)
    else:
        _BLASTN_AUDIT.append(rec)

    if CONFIG["output"].get("blast_progress") and CONFIG["output"].get("log_hits", True):
        ss = "" if same_species is None else (" ✓same" if same_species else " ✗off")
        print(
            f"[HIT {'P' if kind=='protein' else 'N'}] "
            f"db={db} region={region_tag} acc={accession} "
            f"pid≈{rec['percent_identity']}% cov≈{rec['query_coverage']}% "
            f"sp='{hit_species}'{ss} "
            f"def={rec['definition']}"
        )

# =========================
# UniProt helpers
# =========================
def fetch_uniprot_json(acc_or_iso: str) -> Dict:
    """Baixa JSON do UniProt (ACC ou ACC-isoform)."""
    url = f"https://rest.uniprot.org/uniprotkb/{acc_or_iso}.json"
    r = _http_retry("GET", url)
    return r.json()

def fetch_uniprot_fasta(acc: str, include_isoforms: bool=False) -> Dict[str, str]:
    """
    Retorna dict {uniprot_id(ACC ou ACC-isoform): sequence} .
    Se include_isoforms=True, baixa também isoformas.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    if include_isoforms:
        url += "?includeIsoform=true"
    r = _http_retry("GET", url)
    out: Dict[str, str] = {}
    cur_id = None
    buf = []
    for ln in r.text.splitlines():
        if ln.startswith(">"):
            if cur_id and buf:
                out[cur_id] = "".join(buf).strip()
            hdr = ln[1:].strip()
            # headers típicos: sp|P04637|TP53_HUMAN ...  ou sp|P04637-2|TP53_HUMAN Isoform 2 ...
            m = re.search(r"\|([A-Z0-9]+(?:-\d+)?)\|", hdr)
            cur_id = m.group(1) if m else hdr.split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if cur_id and buf:
        out[cur_id] = "".join(buf).strip()
    return out

def uniprot_sequence_from_json(entry: Dict) -> str:
    v = entry.get("sequence", {}).get("value") or ""
    return re.sub(r"\s+", "", v)

def uniprot_taxonomy(entry: Dict) -> Tuple[str, Optional[int]]:
    org = entry.get("organism", {}).get("scientificName") or "?"
    taxid = entry.get("organism", {}).get("taxonId")
    try:
        taxid = int(taxid) if taxid is not None else None
    except Exception:
        taxid = None
    return org, taxid

def extract_features_as_regions(entry: Dict, seq_len: int) -> List[Region]:
    cfg = CONFIG["regions"]
    if not cfg.get("use_uniprot_features", True):
        return []

    incl = set(t.lower() for t in cfg.get("include_feature_types", []))
    point_flank = int(cfg.get("point_flank", 5))   # ±5 por padrão
    min_len = int(cfg.get("default_min_len", 1))
    do_merge = bool(cfg.get("merge_overlaps", True))

    out: List[Region] = []
    for feat in entry.get("features", []) or []:
        ftype = (feat.get("type") or "").lower()
        if incl and ftype not in incl:
            continue

        loc = feat.get("location") or {}
        beg = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if beg is None or end is None:
            continue
        try:
            beg = int(beg); end = int(end)
        except Exception:
            continue

        note = (feat.get("description") or "")[:140].strip()

        # ---- Flanco SOMENTE para feature pontual ----
        if beg == end:
            beg = max(1, beg - point_flank)
            end = min(seq_len, end + point_flank)
        # (Se for intervalo, não aplica flanco)

        if end >= beg and (end - beg + 1) >= min_len:
            out.append(Region(beg, end, ftype.upper(), note))

    # Mescla sobreposições do MESMO tipo, se solicitado
    if do_merge and out:
        out.sort(key=lambda r: (r.tag, r.start_1based, r.end_1based))
        merged: List[Region] = []
        cur = out[0]
        for r in out[1:]:
            if r.tag == cur.tag and r.start_1based <= cur.end_1based + 1:
                # funde intervalos contíguos/overlap do mesmo tipo
                cur.end_1based = max(cur.end_1based, r.end_1based)
            else:
                merged.append(cur)
                cur = r
        merged.append(cur)
        out = merged

    if cfg.get("add_full", True):
        out.insert(0, Region(1, seq_len, "FULL", "FULL"))

    return out

# =========================
# UniProt ⇄ outras bases
# =========================
def uniprot_search_by_refseq_protein(refseq_acc: str) -> Optional[str]:
    """
    Tenta mapear RefSeq Protein → UniProt (ACC canônico).
    Usa REST search por xref:RefSeq:ACC. Retorna ACC (sem isoforma) se achar.
    """
    q_variants = [
        f"xref:RefSeq:{refseq_acc}",
        f"database:RefSeq {refseq_acc}",
    ]
    for q in q_variants:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": q, "format": "json", "size": "1", "fields": "accession"}
        try:
            r = _http_retry("GET", url, params=params)
            j = r.json()
            results = j.get("results", []) or []
            if results:
                return results[0]["primaryAccession"]
        except Exception:
            continue
    return None

# ========== (1) ncbi_gene — choose_best_uniprot_isoform com kwargs inesperado ==========
def choose_best_uniprot_isoform(
    accs: List[str],
    isoform_policy: str = "longest",
    refseq_hint_seq: Optional[str] = None
) -> Optional[str]:
    """
    Escolhe uma isoforma UniProt dentre 'accs'.
    - isoform_policy="longest": pega a mais longa
    - se refseq_hint_seq vier, prioriza ACCs com comprimento igual ao da dica (seq RefSeq)
    """
    if not accs:
        return None

    hint_len = len(refseq_hint_seq) if refseq_hint_seq else None

    def _seq_len_for_acc(acc: str) -> int:
        try:
            fx = fetch_uniprot_fasta(acc, include_isoforms=False)
            # fetch_uniprot_fasta retorna dict {id: seq} — pega a primeira sequência
            if isinstance(fx, dict) and fx:
                return len(next(iter(fx.values())))
            # fallback se alguma implementação retornar str
            if isinstance(fx, str):
                return len(fx)
        except Exception:
            pass
        return 0

    scored: List[Tuple[int, str]] = []
    for acc in accs:
        L = _seq_len_for_acc(acc)
        score = 0
        if hint_len and L == hint_len:
            score += 1_000_000
        if isoform_policy == "longest":
            score += L
        scored.append((score, acc))

    scored.sort(reverse=True)
    return scored[0][1] if scored else accs[0]

# =========================
# PDBe SIFTS (PDB → UniProt)
# =========================
def fetch_pdb_uniprot_mappings(pdb_id: str) -> List[PDBUniProtMap]:
    """
    PDBe SIFTS: PDB -> UniProt mapping
    https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}

    Estrutura (simplificada):
    {
      "<pdb_id>": {
        "UniProt": {
          "<uniprot_acc>": {
            "mappings": [
              {"chain_id": "A", "unp_start": 35, "unp_end": 220, ...},
              ...
            ],
            "sequence_length": 250,
            ...
          },
          ...
        }
      }
    }
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    r = _http_retry("GET", url)
    j = r.json()

    maps: List[PDBUniProtMap] = []
    block = j.get(pdb_id.lower())
    if not block:
        return maps

    uni = block.get("UniProt") or {}
    for acc, obj in uni.items():
        seq_len = None
        try:
            seq_len = int(obj.get("sequence_length") or 0) or None
        except Exception:
            seq_len = None

        for m in (obj.get("mappings") or []):
            # campos robustos, com fallback
            chain = (m.get("chain_id") or m.get("chain") or "?")
            try:
                start = int(m.get("unp_start") or 1)
                end   = int(m.get("unp_end") or start)
            except Exception:
                start, end = 1, 1

            # coverage fornecida pela API ou reconstituída
            cov = m.get("coverage")
            if cov is None:
                if seq_len and end >= start:
                    cov = (end - start + 1) / float(seq_len)
                else:
                    cov = 0.0
            try:
                cov = max(0.0, min(1.0, float(cov)))
            except Exception:
                cov = 0.0

            maps.append(PDBUniProtMap(uniprot_acc=acc, chain=str(chain), coverage=cov))

    return maps

def choose_best_mapping(maps: List[PDBUniProtMap]) -> PDBUniProtMap:
    if not maps: raise ValueError("Sem mapeamentos SIFTS")
    return max(maps, key=lambda m: m.coverage)

# =========================
# Cross-refs UniProt → nuccore/protein
# =========================
def extract_nuccore_and_prot_from_uniprot(entry: Dict) -> Tuple[List[str], List[str]]:
    nuccore: List[str] = []  # NM_, XM_, etc
    prot_ids: List[str] = [] # NP_, XP_, YP_ etc
    for x in entry.get("uniProtKBCrossReferences", []) or []:
        db = (x.get("database") or "").upper()
        xid = (x.get("id") or "").strip()
        if db == "REFSEQ" and xid:
            if re.match(r"^[NXYP]P_[0-9]+\.[0-9]+$", xid):
                prot_ids.append(xid)
        if db in ("REFSEQ","EMBL","GENBANK","DDBJ"):
            for p in x.get("properties", []) or []:
                k = (p.get("key") or "").lower()
                v = (p.get("value") or "").strip()
                if not v: continue
                if "nucleotide" in k or "nucleotid" in k:
                    parts = re.split(r"[;,]\s*", v)
                    for s in parts:
                        if re.match(r"^[A-Z]{1,3}[_\.][A-Za-z0-9_\.]+$", s):
                            nuccore.append(s)
    nuccore = list(dict.fromkeys(nuccore))
    prot_ids = list(dict.fromkeys(prot_ids))
    return nuccore, prot_ids

# =========================
# (2) ncbi_protein — fallbacks RefSeq→UniProt e GenPept features
# =========================
def _strip_refseq_version(acc: str) -> str:
    """
    Remove a versão de um accession RefSeq, se existir.
    Ex.: NP_123.1 -> NP_123 ; XP_054213393.1 -> XP_054213393
    """
    try:
        return acc.split('.')[0]
    except Exception:
        return acc

def map_refseq_protein_to_uniprot(ncbi_prot_acc: str) -> List[str]:
    """
    Mapeia RefSeq Protein (NP_/XP_/YP_/WP_/etc) → UniProt ACC(s) de forma robusta.
    Estratégia (ordem):
      1) xref direto (com versão e sem versão):   xref:RefSeq:ACC
      2) busca por base de dados:                 database:RefSeq AND "ACC"
      3) xref alternativo utilizado em alguns:    xref:RefSeq_Protein:ACC
      4) Fallback por busca livre (texto):        "ACC" e ACC (com/sem versão)
         - cobre casos TrEMBL sem xref curado, mas onde o accession aparece no texto/descrição
    - Nunca levanta exceção por 400/5xx: ignora a tentativa e prossegue.
    - Retorna lista deduplicada preservando ordem.
    - Loga consultas tentadas quando não houver resultado e blast_progress=True.
    """
    base = _strip_refseq_version(ncbi_prot_acc)
    tried_queries: List[str] = []

    def _try_query(q: str) -> List[str]:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": q, "fields": "accession", "format": "json", "size": "50"}
        tried_queries.append(f"{url}?query={q}")
        try:
            r = _http_retry("GET", url, params=params)
            j = r.json() or {}
            out: List[str] = []
            for it in (j.get("results") or []):
                acc = it.get("primaryAccession") or it.get("uniProtkbId")
                if acc:
                    out.append(str(acc))
            return out
        except Exception:
            # qualquer erro (400/5xx/network) -> retorna vazio e segue
            return []

    results: List[str] = []

    # ---- (1) xref direto (com e sem versão) ----
    results += _try_query(f"xref:RefSeq:{ncbi_prot_acc}")
    if not results:
        results += _try_query(f"xref:RefSeq:{base}")

    # ---- (2) database:RefSeq AND "<acc>" (com e sem versão) ----
    if not results:
        results += _try_query(f'database:RefSeq AND "{ncbi_prot_acc}"')
    if not results:
        results += _try_query(f'database:RefSeq AND "{base}"')

    # ---- (3) xref alternativo RefSeq_Protein (com e sem versão) ----
    if not results:
        results += _try_query(f"xref:RefSeq_Protein:{ncbi_prot_acc}")
    if not results:
        results += _try_query(f"xref:RefSeq_Protein:{base}")

    # ---- (4) Fallback por busca livre (texto) ----
    # primeiro com aspas (match mais estrito), depois sem aspas
    if not results:
        results += _try_query(f'"{ncbi_prot_acc}"')
    if not results:
        results += _try_query(f'"{base}"')
    if not results:
        results += _try_query(ncbi_prot_acc)
    if not results:
        results += _try_query(base)

    # Dedup preservando ordem
    results = list(dict.fromkeys(results))

    if not results and CONFIG.get("output", {}).get("blast_progress", False):
        print(f"[DIAG] UniProt search sem resultados para {ncbi_prot_acc}. Consultas tentadas: {tried_queries}")

    return results

def extract_regions_from_genpept_features(refseq_prot_acc: str) -> List[Region]:
    """
    Fallback de features a partir do GenPept (db=protein, rettype=gp).
    Converte as features para Region (somente intervalares; sítios pontuais recebem flanco padrão de 5 aa).
    """
    # baixa o registro GenPept
    handle = Entrez.efetch(db="protein", id=refseq_prot_acc, rettype="gp", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    seq_len = len(record.seq) if record and record.seq else 0
    out: List[Region] = []
    flank = 5

    # mapeamento simples de tipos (GenPept → tags)
    interval_types = {
        "Region": "REGION",
        "Site": "SITE",
        "Domain": "DOMAIN",
        "Repeat": "REPEAT",
        "Motif": "MOTIF",
        "Transmembrane": "TRANSMEMBRANE",
        "Intramembrane": "INTRAMEMBRANE",
        "Coiled-coil": "COILED COIL",
        "Zinc finger": "ZINC FINGER",
        "Signal peptide": "Signal peptide".upper(),
        "Transit peptide": "TRANSIT PEPTIDE",
        "Propeptide": "PROPEPTIDE",
        "peptide": "PEPTIDE",
        "Chain": "CHAIN",
        "Topological domain": "TOPOLOGICAL DOMAIN",
        "Active site": "ACTIVE SITE",
        "Metal binding": "METAL BINDING",
        "Binding site": "BINDING SITE",
        "Glycosylation": "GLYCOSYLATION SITE",
        "Lipidation": "LIPIDATION",
        "Modified residue": "MODIFIED RESIDUE",
        "Disulfide bond": "DISULFIDE BOND",
        "Cross-link": "CROSS-LINK",
    }

    for feat in record.features:
        if feat.type not in ("Region", "Site", "Domain", "Repeat", "Motif",
                             "Transmembrane", "Intramembrane", "Coiled-coil",
                             "Zinc finger", "Signal peptide", "Transit peptide",
                             "Propeptide", "peptide", "Chain", "Topological domain",
                             "Active site", "Metal binding", "Binding site",
                             "Glycosylation", "Lipidation", "Modified residue",
                             "Disulfide bond", "Cross-link"):
            continue

        tag = interval_types.get(feat.type, feat.type.upper())
        note = ""
        try:
            note = (feat.qualifiers.get("note", [""])[0] or "")[:140]
        except Exception:
            pass

        # localização
        try:
            beg = int(feat.location.nofuzzy_start) + 1   # 1-based
            end = int(feat.location.nofuzzy_end)
        except Exception:
            continue

        if beg == end:
            # flanco apenas para sítios pontuais
            beg = max(1, beg - flank)
            end = min(seq_len, end + flank)

        if end >= beg and (end - beg + 1) >= int(CONFIG["regions"]["default_min_len"]):
            out.append(Region(beg, end, tag, note))

    # FULL sempre presente
    out.insert(0, Region(1, seq_len, "FULL", "FULL"))
    return out

# =========================
# NCBI: protein→nuccore (elink) e GenBank → CDS correta
# =========================
def elink_protein_to_nuccore_accs(prot_acc: str) -> List[str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {"dbfrom":"protein","db":"nuccore","id":prot_acc,"retmode":"json"}
    if NCBI_API_KEY: params["api_key"] = NCBI_API_KEY
    if NCBI_EMAIL:   params["email"] = NCBI_EMAIL
    r = _http_retry("GET", url, params=params)
    j = r.json()
    out = []
    try:
        linksets = j["linksets"][0].get("linksetdbs", []) or []
        ids = []
        for ls in linksets:
            ids.extend(ls.get("links", []) or [])
        if ids:
            url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            p2 = {"db":"nuccore","id":",".join(str(x) for x in ids),"rettype":"acc","retmode":"text"}
            if NCBI_API_KEY: p2["api_key"] = NCBI_API_KEY
            if NCBI_EMAIL:   p2["email"] = NCBI_EMAIL
            r2 = _http_retry("GET", url2, params=p2)
            for line in r2.text.splitlines():
                acc = line.strip()
                if acc: out.append(acc)
    except Exception:
        pass
    return list(dict.fromkeys(out))

def pick_cds_from_genbank(record, prot_seq: str, candidate_protein_ids: List[str]) -> Tuple[Optional[str], Optional[str], Optional[Dict]]:
    """Retorna (cds_nt, translation, feat_meta) da CDS que:
       - tem /protein_id em candidate_protein_ids OU
       - tem /translation == prot_seq (exata, ignorando '*' final)
    """
    best = None
    for feat in record.features:
        if feat.type != "CDS": continue
        q = feat.qualifiers or {}
        protein_id = (q.get("protein_id",[None])[0]) or None
        translation = (q.get("translation",[None])[0]) or None
        codon_start = int((q.get("codon_start",[1])[0]) or 1)
        transl_table = int((q.get("transl_table",[1])[0]) or 1)

        spliced = str(feat.extract(record.seq)).upper()
        if codon_start in (2,3):
            spliced = spliced[codon_start-1:]
        if len(spliced) % 3 != 0:
            spliced = spliced[:len(spliced)//3*3]

        try:
            tr = str(Seq(spliced).translate(table=transl_table, to_stop=False))
        except Exception:
            tr = None
        tr_nostop = tr[:-1] if (tr and tr.endswith("*")) else tr

        match_protid = (protein_id in set(candidate_protein_ids)) if protein_id else False
        prot_seq_nostop = prot_seq[:-1] if prot_seq.endswith("*") else prot_seq
        match_translation = (translation == prot_seq) or (tr_nostop == prot_seq_nostop)

        score = (1 if match_protid else 0, 1 if match_translation else 0)
        if score == (0,0):
            continue
        meta = {"protein_id": protein_id, "codon_start": codon_start, "transl_table": transl_table}
        cand = (spliced, translation, meta, score)
        if best is None or cand[3] > best[3]:
            best = cand
    if not best:
        return None, None, None
    return best[0], best[1], best[2]

def find_cds_strict(nuccore_accs: List[str], prot_seq: str, candidate_protein_ids: List[str]) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[Dict]]:
    """Percorre nuccores e retorna a primeira CDS consistente.
       Out: (picked_nuccore, cds_nt, translation, feat_meta)
    """
    for acc in nuccore_accs:
        try:
            if CONFIG["output"]["blast_progress"]:
                print(f"[DIAG] Verificando CDS em {acc} …")
            handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            cds_nt, translation, meta = pick_cds_from_genbank(record, prot_seq, candidate_protein_ids)
            if not cds_nt:
                continue
            prot_seq_nostop = prot_seq[:-1] if prot_seq.endswith("*") else prot_seq
            if len(cds_nt) != 3*len(prot_seq_nostop):
                tr = str(Seq(cds_nt).translate(to_stop=False))
                if tr.endswith("*"):
                    tr = tr[:-1]
                if tr != prot_seq_nostop:
                    continue
            return acc, cds_nt, translation, meta
        except Exception as e:
            if CONFIG["output"]["blast_progress"]:
                print(f"[DIAG]  - falha em {acc}: {e}")
            continue
    return None, None, None, None

# =========================
# TaxID helpers — mesma espécie p/ AA e NT
# =========================
def _get_tax_from_genpept(prot_acc: str) -> Tuple[Optional[int], Optional[str]]:
    """
    Busca taxid e nome do organismo direto do GenPept (RefSeq Protein).
    - taxid: em feature 'source' (qualifier 'db_xref': 'taxon:<id>')
    - organismo: 'organism' nas annotations ou qualifier 'organism' em 'source'
    """
    tax = None
    org = None
    try:
        handle = Entrez.efetch(db="protein", id=prot_acc, rettype="gp", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception:
        return None, None

    try:
        org = record.annotations.get("organism")
    except Exception:
        org = None
    try:
        for feat in record.features:
            if feat.type == "source":
                for x in feat.qualifiers.get("db_xref", []):
                    if x.startswith("taxon:"):
                        tax = int(x.split(":")[1])
                if not org:
                    org = (feat.qualifiers.get("organism", [""])[0] or None)
    except Exception:
        pass
    return tax, org

def ensure_blast_same_species_filters_from_taxid(taxid: Optional[int], species: Optional[str]=None):
    """
    Se same_species_only=True e houver taxid, setamos entrez_query para AA e NT.
    Também fixa o alvo esperado p/ auditoria (espécie).
    """
    global _EXPECTED_TAXID, _EXPECTED_SPECIES
    if species:
        _EXPECTED_SPECIES = species
    if taxid and CONFIG.get("blast", {}).get("same_species_only", False):
        q = f"txid{int(taxid)}[ORGN]"
        CONFIG["blast"]["protein"]["entrez_query"] = q
        CONFIG["blast"]["nt"]["entrez_query"] = q
        _EXPECTED_TAXID = int(taxid)
        if CONFIG.get("output", {}).get("blast_progress", False):
            print(f"[BLAST] Filtro por espécie aplicado: {q}")

# =========================
# BLAST por região
# =========================
def blastp_region(seq_window: str, entrez_query: Optional[str], params: Dict) -> Dict[int, set]:
    if not seq_window or not re.search(r"[A-Z]", seq_window):
        return {}
    handle = NCBIWWW.qblast(
        program="blastp",
        database=params.get("db", "nr"),
        sequence=seq_window,
        hitlist_size=int(params.get("hitlist_size", 25)),
        expect=float(params.get("expect", 1e-5)),
        entrez_query=entrez_query,
        service="plain",
        megablast=False
    )
    rec = NCBIXML.read(handle)
    handle.close()
    variants: Dict[int, set] = {}
    L = len(seq_window)
    min_id = float(params.get("min_identity", 0.97))
    min_cov = float(params.get("min_query_coverage", 0.95))
    region_tag = params.get("audit_region_tag", "?")
    db_used = params.get("db", "nr")

    for aln in rec.alignments:
        for hsp in aln.hsps:
            qseq = hsp.query
            hseq = hsp.sbjct
            q_aligned = sum(1 for c in qseq if c != '-')
            if q_aligned == 0:
                continue
            cov = q_aligned / L
            if cov < min_cov:
                continue
            ident = hsp.identities / q_aligned
            if ident < min_id:
                continue

            # ===== Auditoria: log/CSV do hit aceito =====
            _audit_hit(
                kind="protein",
                db=db_used,
                region_tag=region_tag,
                accession=getattr(aln, "accession", "") or "",
                hit_def=getattr(aln, "hit_def", "") or "",
                pid=ident,
                cov=cov,
                hsp_query_len=q_aligned
            )

            # ===== variantes / colchetes =====
            qpos = 0
            for qc, hc in zip(qseq, hseq):
                if qc != '-':
                    qpos += 1
                    if hc != '-' and qc != hc:
                        variants.setdefault(qpos, set()).add(hc)
    return variants

def blastp_region_multi(seq_window: str, entrez_query: Optional[str], params: Dict) -> Dict[int, set]:
    """
    Executa BLASTp em múltiplos bancos (params['dbs'] ou [params['db']]) e
    retorna a UNIÃO dos conjuntos de alternativas por posição.
    Respeita min_identity / min_query_coverage do próprio 'params'.
    """
    dbs = params.get("dbs")
    if not dbs:
        dbs = [params.get("db", "nr")]
    union_variants: Dict[int, set] = {}

    for db in dbs:
        if CONFIG["output"]["blast_progress"]:
            print(f"[BLASTp] Rodando em {db}… (window={len(seq_window)} aa)")
        local = blastp_region(seq_window, entrez_query, {**params, "db": db})
        for pos, alts in local.items():
            union_variants.setdefault(pos, set()).update(alts)

    return union_variants

def blastn_region(nt_window: str, db: str, params: Dict, is_mrna: bool=False) -> Dict[int, set]:
    """
    BLASTn por janela. Trabalhamos em DNA (T). Para mRNA, convertemos U->T antes do BLAST.
    Agora suporta entrez_query (mesmo filtro de espécie do AA).
    Se CONFIG["blast"]["log_offspecies_preview"] = True e houver entrez_query,
    roda uma 2ª passada só para LOG (sem afetar variantes) sem filtro de espécie.
    """
    if not nt_window or not re.search(r"[ACGTU]", nt_window, re.I):
        return {}

    query_dna = nt_window.upper().replace("U", "T")
    entrez_q = params.get("entrez_query")
    if CONFIG["output"]["blast_progress"]:
        print(f"[BLASTn] Rodando em {db} (len={len(query_dna)} nt) "
              f"{'(mesma espécie)' if entrez_q else '(todas as espécies)'}…")

    handle = NCBIWWW.qblast(
        program="blastn",
        database=db,
        sequence=query_dna,
        hitlist_size=int(params.get("hitlist_size", 25)),
        expect=float(params.get("expect", 1e-10)),
        service="plain",
        megablast=bool(params.get("megablast", True)),
        entrez_query=entrez_q
    )
    rec = NCBIXML.read(handle); handle.close()

    variants: Dict[int, set] = {}
    L = len(query_dna)
    min_id = float(params.get("min_identity", 0.98))
    min_cov = float(params.get("min_query_coverage", 0.95))
    region_tag = params.get("audit_region_tag", "?")

    def _log_hits(rec_obj, level_label: str):
        for aln in rec_obj.alignments:
            # melhor HSP por alinhamento
            best_cov = 0.0; best_pid = 0.0
            best_hsp = None
            for hsp in aln.hsps:
                q_aligned = sum(1 for c in hsp.query if c != '-')
                cov = (q_aligned / L) if L else 0.0
                pid = (hsp.identities / q_aligned) if q_aligned else 0.0
                if (pid, cov) > (best_pid, best_cov):
                    best_pid, best_cov, best_hsp = pid, cov, hsp
            _audit_hit(
                kind="nucleotide",
                db=f"{db}{level_label}",
                region_tag=region_tag,
                accession=getattr(aln, "accession", "") or "",
                hit_def=getattr(aln, "hit_def", "") or "",
                pid=best_pid,
                cov=best_cov,
                hsp_query_len=best_hsp.align_length if best_hsp else 0
            )

    # ---- HITS usados (com ou sem filtro de espécie) ----
    if CONFIG["output"]["blast_progress"]:
        _log_hits(rec, "")

    for aln in rec.alignments:
        for hsp in aln.hsps:
            qseq = hsp.query.upper()
            hseq = hsp.sbjct.upper()
            q_aligned = sum(1 for c in qseq if c != '-')
            if not q_aligned or (q_aligned / L) < min_cov:
                continue
            ident = (hsp.identities / q_aligned)
            if ident < min_id:
                continue
            qpos = 0
            for qc, hc in zip(qseq, hseq):
                if qc != '-':
                    qpos += 1
                    if hc != '-' and qc != hc:
                        variants.setdefault(qpos, set()).add(hc)

    # ---- PREVIEW off-species (apenas log) ----
    try:
        if CONFIG["blast"].get("log_offspecies_preview", False) and entrez_q:
            if CONFIG["output"]["blast_progress"]:
                print("[BLASTn] Preview off-species (somente log; não altera variantes)…")
            handle2 = NCBIWWW.qblast(
                program="blastn",
                database=db,
                sequence=query_dna,
                hitlist_size=min(int(params.get("hitlist_size", 25)), 5),
                expect=float(params.get("expect", 1e-10)),
                service="plain",
                megablast=bool(params.get("megablast", True)),
                entrez_query=None
            )
            rec2 = NCBIXML.read(handle2); handle2.close()
            _log_hits(rec2, " PREVIEW")
    except Exception as e:
        if CONFIG["output"]["blast_progress"]:
            print(f"[WARN] Preview off-species falhou: {e}")

    return variants

# =========================
# Effatha builders
# =========================
def effatha_aa(seq: str, region: Region, local_variants: Dict[int, set]) -> str:
    s = region.start_1based
    e = region.end_1based
    buf = []
    for i, pos in enumerate(range(s, e+1), start=1):
        ref = seq[pos-1]
        alts = sorted(x for x in local_variants.get(i, set()) if x != ref)
        if alts:
            buf.append("[" + "/".join([ref] + alts) + "]")
        else:
            buf.append(ref)
    return "".join(buf)

def effatha_nt(nt_seq: str, nt_variants: Dict[int, set], use_u: bool=False) -> str:
    buf = []
    for i, ref in enumerate(nt_seq, start=1):
        ref_char = ref
        alts = sorted(x for x in nt_variants.get(i, set()) if x != ref_char)
        if use_u:
            ref_char = 'U' if ref_char.upper() == 'T' else ref_char
            alts = [('U' if a.upper() == 'T' else a) for a in alts]
        if alts:
            buf.append("[" + "/".join([ref_char] + alts) + "]")
        else:
            buf.append(ref_char)
    return "".join(buf)

# =========================
# Núcleo: normalizar entrada → UniProt → features → BLAST → NT
# =========================
def core_pipeline_using_uniprot(target_uniprot: str,
                                refseq_hint_seq: Optional[str]=None,
                                mapping_info: Optional[Dict]=None,
                                extra_prot_ids: Optional[List[str]]=None):
    """
    target_uniprot pode ser ACC (P04637) ou isoforma (P04637-2).
    refseq_hint_seq: sequência de proteína (ex.: da entrada RefSeq) usada para tentar escolher isoforma.
    extra_prot_ids: lista de protein_ids RefSeq (para ajudar a travar a CDS correta).
    """
    # Logs de início (4) uniprot — “não roda/não loga”
    print(f"[INFO] UniProt {target_uniprot}: iniciando coleta de sequência, features e CDS…")

    entry = fetch_uniprot_json(target_uniprot)
    seq = uniprot_sequence_from_json(entry)
    org, taxid = uniprot_taxonomy(entry)

    # Filtro de mesma espécie (AA e NT) + alvo esperado p/ auditoria
    ensure_blast_same_species_filters_from_taxid(taxid, species=org)

    regs = extract_features_as_regions(entry, len(seq))

    # ===== NOVO: fallback automático para GenPept quando UniProt vier sem features reais =====
    if CONFIG["regions"].get("fallback_genpept_if_uniprot_featureless", True):
        only_full_or_empty = (not regs) or all(r.tag == "FULL" for r in regs)
        if only_full_or_empty:
            # tentar obter protein_ids RefSeq do próprio UniProt + extras do caller
            _nuccs_from_up, prot_ids_from_up = extract_nuccore_and_prot_from_uniprot(entry)
            prot_candidates: List[str] = list(dict.fromkeys((extra_prot_ids or []) + prot_ids_from_up))
            if prot_candidates:
                print(f"[INFO] UniProt sem features (ou somente FULL). Tentando fallback GenPept usando {prot_candidates[0]} …")
                try:
                    regs_gp = extract_regions_from_genpept_features(prot_candidates[0])
                    if regs_gp:
                        regs = regs_gp
                        print("[OK] Features do GenPept aplicadas como fallback.")
                except Exception as e:
                    print(f"[WARN] Fallback GenPept falhou: {e}")

    if not regs:
        regs = [Region(1, len(seq), "FULL", "FULL")]

    # ===== AA por região (BLASTp → Effatha) =====
    aa_by_region: List[Tuple[str,int,int,str,str]] = []
    blast_positions_count = 0
    if CONFIG["blast"]["enable"]:
        bp = CONFIG["blast"]["protein"]
        entrez_q = bp.get("entrez_query")
        for r in regs:
            frag = seq[r.start_1based-1:r.end_1based]
            if CONFIG["output"]["blast_progress"]:
                dbs_log = bp.get("dbs") or [bp.get("db","nr")]
                print(f"[BLASTp] {r.tag} ({len(frag)} aa) em {dbs_log} (mesmo taxon via entrez_query={entrez_q})…")
            try:
                bp_with_tag = {**bp, "audit_region_tag": r.tag}
                local_vars = blastp_region_multi(frag, entrez_q, bp_with_tag)
            except Exception as e:
                if CONFIG["output"]["blast_progress"]:
                    print(f"[WARN] BLASTP falhou em {r.tag} {r.start_1based}-{r.end_1based}: {e}")
                local_vars = {}
            blast_positions_count += len(local_vars)
            aa_eff = effatha_aa(seq, r, local_vars)
            aa_by_region.append((r.tag, r.start_1based, r.end_1based, aa_eff, r.note))
    else:
        for r in regs:
            frag = seq[r.start_1based-1:r.end_1based]
            aa_by_region.append((r.tag, r.start_1based, r.end_1based, frag, r.note))

    # ===== NT por região (CDS real → fatiamento → BLASTn → Effatha) =====
    nt_segments: List[Tuple[str,int,int,str,str,str]] = []
    picked_nuccore = None
    if CONFIG["cds_mapping"]["enable"]:
        nuccore_cands, prot_ids_from_up = extract_nuccore_and_prot_from_uniprot(entry)
        extra_from_prot = []
        for pid in (extra_prot_ids or []):
            extra_from_prot += elink_protein_to_nuccore_accs(pid)
        for pid in (prot_ids_from_up or []):
            extra_from_prot += elink_protein_to_nuccore_accs(pid)
        nuccore_all = list(dict.fromkeys(nuccore_cands + extra_from_prot))
        if CONFIG["output"]["blast_progress"]:
            print(f"[DIAG] nuccore candidatos: {nuccore_all or '—'} ; protein_ids: {(extra_prot_ids or []) + prot_ids_from_up or ['—']}")

        picked_nuccore, cds_nt, _translation, _meta = find_cds_strict(nuccore_all, seq, (extra_prot_ids or []) + prot_ids_from_up)
        if picked_nuccore and cds_nt:
            bn = CONFIG["blast"]["nt"]
            for r in regs:
                dna = cds_nt[(r.start_1based-1)*3 : r.end_1based*3]
                mrna = dna.replace("T","U")
                dna_vars = {}
                mrna_vars = {}
                if CONFIG["blast"]["enable"]:
                    try:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[BLASTn] DNA {r.tag} ({len(dna)} nt) em {bn.get('dna_db','nt')}…")
                        bn_with_tag = {**bn, "audit_region_tag": r.tag}
                        dna_vars = blastn_region(dna, bn.get("dna_db","nt"), bn_with_tag, is_mrna=False)
                    except Exception as e:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[WARN] BLASTN DNA falhou em {r.tag} {r.start_1based}-{r.end_1based}: {e}")
                    try:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[BLASTn] mRNA {r.tag} ({len(mrna)} nt) em {bn.get('rna_db','refseq_rna')}…")
                        bn_with_tag = {**bn, "audit_region_tag": r.tag}
                        mrna_vars = blastn_region(mrna, bn.get("rna_db","refseq_rna"), bn_with_tag, is_mrna=True)
                    except Exception as e:
                        if CONFIG["output"]["blast_progress"]:
                            print(f"[WARN] BLASTN mRNA falhou em {r.tag} {r.start_1based}-{r.end_1based}: {e}")

                dna_eff  = effatha_nt(dna, dna_vars, use_u=False)
                mrna_eff = effatha_nt(mrna, mrna_vars, use_u=True)
                nt_segments.append((r.tag, r.start_1based, r.end_1based, dna_eff, mrna_eff, r.note))
        else:
            if CONFIG["output"]["blast_progress"]:
                print("[DIAG] CDS real não mapeada — NT ficará ausente para estas regiões.")

    # ===== Contexto e artifacts =====
    protein_meta = {"accession": target_uniprot, "length": len(seq), "organism": org, "taxid": taxid}
    if mapping_info:
        protein_meta["mapping"] = mapping_info

    region_cards = []
    for r in regs:
        region_cards.append({
            "tag": r.tag, "start": r.start_1based, "end": r.end_1based, "length": r.length, "note": r.note,
            "ligand_name": None, "ligand_chebi": None, "ligand_is_metal": None,
            "pct_var_obs": 0.0, "pct_var_fil": 0.0, "max_af": None, "clin_counts": {},
            "pos_with_manual_ev": 0, "sift_tol": 0, "sift_del": 0, "poly_benign": 0, "poly_damaging": 0,
            "has_deletion": False,
        })

    consolidated = {
        "regions": [{"tag": r.tag, "start": r.start_1based, "end": r.end_1based} for r in regs],
        "region_cards": region_cards,
        "aa_regions": [{"tag": t, "start": s, "end": e, "aa": a} for (t,s,e,a,_) in aa_by_region],
        "aa_regions_obs": [],
        "aa_regions_filtered": [],
        "nt_segments": [
            {"tag": t, "start": s, "end": e, "dna": dna, "mrna": mrna, "note": note}
            for (t,s,e,dna,mrna,note) in nt_segments
        ],
    }

    return seq, protein_meta, consolidated

# =========================
# Entradas
# =========================
def run_from_uniprot(uniprot_acc: str):
    seq, meta, consolidated = core_pipeline_using_uniprot(uniprot_acc)
    context = {
        "protein": meta,
        "layers": {
            "uniprot_variant_enabled": True,
            "proteins_variation_enabled": True,
            "blast_enabled": bool(CONFIG["blast"]["enable"]),
        },
        "blast": {
            "variant_positions_blast": 0,  # contador fino pode ser adicionado se necessário
            "variant_positions_uniprot": 0,
            "variant_positions_proteins": 0,
            "variant_positions_filtered": 0,
        },
        **consolidated,
        # Abas: tudo nasce da visão UniProt; NT (CDS) também fica disponível na aba gene
        "sources": {
            "uniprot": consolidated,
            "pdb": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_protein": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": consolidated["nt_segments"]},
        }
    }
    _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)

def run_from_pdb(pdb_id: str):
    pdb = (pdb_id or "").strip()
    if not pdb:
        print("[ERRO] PDB vazio.")
        return
    print(f"Descobrindo UniProt para PDB {pdb} via SIFTS…")
    maps = fetch_pdb_uniprot_mappings(pdb)
    if not maps:
        print("[ERRO] Nenhum mapeamento UniProt encontrado via SIFTS.")
        return
    best = choose_best_mapping(maps)
    chosen = best.uniprot_acc
    print(f"Selecionado UniProt {chosen} (cadeia {best.chain}, cobertura≈{best.coverage:.1%})")

    # setar filtro por taxid (quando disponível) — aplica em AA e NT + auditoria alvo
    org_tax_entry = fetch_uniprot_json(chosen)
    _org, _taxid = uniprot_taxonomy(org_tax_entry)
    ensure_blast_same_species_filters_from_taxid(_taxid, species=_org)

    mapping_info = {"input_source":"pdb","input_id":pdb_id,"mapped_uniprot":chosen,"notes":"SIFTS (PDBe)"}
    seq, meta, consolidated = core_pipeline_using_uniprot(chosen, mapping_info=mapping_info)
    context = {
        "protein": meta,
        "layers": {
            "uniprot_variant_enabled": True,
            "proteins_variation_enabled": True,
            "blast_enabled": bool(CONFIG["blast"]["enable"]),
        },
        "blast": {"variant_positions_blast": 0, "variant_positions_uniprot": 0, "variant_positions_proteins": 0, "variant_positions_filtered": 0},
        **consolidated,
        "sources": {
            "uniprot": consolidated,
            "pdb": consolidated,  # exibir igual na aba PDB (mesma sequência/regions), anotando mapping em meta
            "ncbi_protein": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
            "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": consolidated["nt_segments"]},
        }
    }
    _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)

def run_from_ncbi_protein(refseq_acc: str):
    acc = (refseq_acc or "").strip()
    if not acc:
        raise ValueError("RefSeq Protein ACC vazio.")

    # Sequência da proteína (RefSeq)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db":"protein","id":acc,"rettype":"fasta","retmode":"text"}
    if NCBI_API_KEY: params["api_key"] = NCBI_API_KEY
    if NCBI_EMAIL:   params["email"] = NCBI_EMAIL
    r = _http_retry("GET", url, params=params)
    refseq_seq = "".join([ln.strip() for ln in r.text.splitlines() if not ln.startswith(">")])
    if not refseq_seq:
        raise RuntimeError("Falha ao obter sequência proteica do RefSeq.")

    # NEW: Descobrir taxid/espécie pelo próprio GenPept (mesma espécie em AA e NT mesmo sem UniProt)
    gp_tax, gp_org = _get_tax_from_genpept(acc)
    ensure_blast_same_species_filters_from_taxid(gp_tax, species=gp_org)

    # (2A) Mapear RefSeq→UniProt via busca no UniProt (reversa)
    up_accs = map_refseq_protein_to_uniprot(acc)

    if up_accs:
        # escolher melhor isoforma/ACC conforme política e dica de sequência
        up = choose_best_uniprot_isoform(up_accs, isoform_policy="longest", refseq_hint_seq=refseq_seq)
        # Rodar pipeline UniProt, informando o protein_id RefSeq para travar CDS
        seq, meta, consolidated = core_pipeline_using_uniprot(
            up,
            refseq_hint_seq=refseq_seq,
            mapping_info={"input_source":"ncbi_protein","input_id":acc,"mapped_uniprot":up},
            extra_prot_ids=[acc]
        )

        context = {
            "protein": meta,
            "layers": {"uniprot_variant_enabled": True, "proteins_variation_enabled": True, "blast_enabled": bool(CONFIG["blast"]["enable"])},
            "blast": {"variant_positions_blast": 0, "variant_positions_uniprot": 0, "variant_positions_proteins": 0, "variant_positions_filtered": 0},
            **consolidated,
            "sources": {
                "uniprot": consolidated,                 # regiões funcionais nas abas UniProt e NCBI Protein
                "pdb": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": []},
                "ncbi_protein": consolidated,
                "ncbi_gene": {"regions": [], "region_cards": [], "aa_regions": [], "aa_regions_obs": [], "aa_regions_filtered": [], "nt_segments": consolidated["nt_segments"]},
            }
        }
        _write_artifacts_and_context(seq, consolidated["aa_regions"], consolidated["nt_segments"], context)
        return

    # (2B) Fallback: extrair features do GenPept quando não há mapeamento UniProt
    print("[AVISO] Não foi possível mapear RefSeq→UniProt pelo xref. Usando features do GenPept.")
    regs = extract_regions_from_genpept_features(acc)
    seq = refseq_seq

    # BLAST AA com logs de progresso
    aa_by_region = []
    if CONFIG["blast"]["enable"]:
        bp = CONFIG["blast"]["protein"]
        for rgn in regs:
            frag = seq[rgn.start_1based-1:rgn.end_1based]
            if CONFIG["output"]["blast_progress"]:
                dbs_log = bp.get("dbs") or [bp.get("db","nr")]
                print(f"[BLASTp] {rgn.tag} ({len(frag)} aa) em {dbs_log} (mesmo taxon via entrez_query={bp.get('entrez_query')})…")
            try:
                bp_with_tag = {**bp, "audit_region_tag": rgn.tag}
                local_vars = blastp_region_multi(frag, bp.get("entrez_query"), bp_with_tag)
            except Exception as e:
                if CONFIG["output"]["blast_progress"]:
                    print(f"[WARN] BLASTP falhou: {e}")
                local_vars = {}
            aa_by_region.append({"tag":rgn.tag,"start":rgn.start_1based,"end":rgn.end_1based,"aa": effatha_aa(seq, rgn, local_vars)})
    else:
        for rgn in regs:
            frag = seq[rgn.start_1based-1:rgn.end_1based]
            aa_by_region.append({"tag":rgn.tag,"start":rgn.start_1based,"end":rgn.end_1based,"aa": frag})

    # NT a partir da própria RefSeq protein→nuccore
    nt_segments = []
    if CONFIG["cds_mapping"]["enable"]:
        picked, cds_nt, _t, _m = find_cds_strict(elink_protein_to_nuccore_accs(acc), seq, [acc])
        if picked and cds_nt:
            bn = CONFIG["blast"]["nt"]
            for rgn in regs:
                dna = cds_nt[(rgn.start_1based-1)*3 : rgn.end_1based*3]
                mrna = dna.replace("T","U")
                try:
                    if CONFIG["blast"]["enable"] and CONFIG["output"]["blast_progress"]:
                        print(f"[BLASTn] DNA {rgn.tag} ({len(dna)} nt) em {bn.get('dna_db','nt')}…")
                    bn_with_tag = {**bn, "audit_region_tag": rgn.tag}
                    dna_vars = blastn_region(dna, bn.get("dna_db","nt"), bn_with_tag, is_mrna=False) if CONFIG["blast"]["enable"] else {}
                except Exception as e:
                    dna_vars = {}
                    print(f"[WARN] BLASTN DNA falhou: {e}")
                try:
                    if CONFIG["blast"]["enable"] and CONFIG["output"]["blast_progress"]:
                        print(f"[BLASTn] mRNA {rgn.tag} ({len(mrna)} nt) em {bn.get('rna_db','refseq_rna')}…")
                    bn_with_tag = {**bn, "audit_region_tag": rgn.tag}
                    mrna_vars = blastn_region(mrna, bn.get("rna_db","refseq_rna"), bn_with_tag, is_mrna=True) if CONFIG["blast"]["enable"] else {}
                except Exception as e:
                    mrna_vars = {}
                    print(f"[WARN] BLASTN mRNA falhou: {e}")
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
        "aa_regions": aa_by_region,
        "aa_regions_obs": [],
        "aa_regions_filtered": [],
        "nt_segments": nt_segments,
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

    # Buscar RefSeq Proteins do gene
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    p = {"db":"protein","term":q,"retmax":"100","sort":"slen","retmode":"json"}
    if NCBI_API_KEY: p["api_key"] = NCBI_API_KEY
    if NCBI_EMAIL:   p["email"] = NCBI_EMAIL
    r = _http_retry("GET", url, params=p)
    ids = (r.json().get("esearchresult",{}).get("idlist") or [])
    if not ids: raise RuntimeError("Nenhum RefSeq Protein encontrado.")

    url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    p2 = {"db":"protein","id":",".join(ids),"rettype":"acc","retmode":"text"}
    if NCBI_API_KEY: p2["api_key"] = NCBI_API_KEY
    if NCBI_EMAIL:   p2["email"] = NCBI_EMAIL
    r2 = _http_retry("GET", url2, params=p2)
    accs = [ln.strip() for ln in r2.text.splitlines() if ln.strip()]
    if not accs: raise RuntimeError("Falha ao obter ACCs de protein.")

    # Escolher a maior isoforma (ou política configurável)
    refseq_acc = max(accs, key=lambda a: len(a))
    return run_from_ncbi_protein(refseq_acc)

# =========================
# Saídas
# =========================
def _write_artifacts_and_context(
    prot_seq: str,
    aa_regions: List[Dict[str, Any]],
    nt_segments: List[Dict[str, Any]],
    context: Dict[str, Any],
):
    # report
    rp = CONFIG["output"]["report_txt"]
    os.makedirs(os.path.dirname(rp), exist_ok=True)
    with io.open(rp, "w", encoding="utf-8") as fh:
        fh.write("# Effatha — Gene & Protein Analyzer (regiões funcionais)\n\n")
        pmeta = context.get("protein", {})
        fh.write(f"Proteína: {pmeta.get('accession','?')}\n")
        fh.write(f"Organismo: {pmeta.get('organism','?')} (taxid={pmeta.get('taxid','?')})\n")
        fh.write(f"Tamanho (AA): {len(prot_seq)} aa\n\n")

        fh.write("## Aminoácidos por região (sintaxe Effatha)\n")
        for a in aa_regions:
            fh.write(f"\n### {a['tag']} {a['start']}-{a['end']}\n")
            fh.write(a["aa"] + "\n")

        fh.write("\n\n## Segmentos NT (DNA/mRNA; sintaxe Effatha por base)\n")
        if nt_segments:
            for s in nt_segments:
                fh.write(f"\n### {s['tag']} {s['start']}-{s['end']}\n")
                fh.write("DNA:\n")
                fh.write(s["dna"] + "\n")
                fh.write("mRNA:\n")
                fh.write(s["mrna"] + "\n")
        else:
            fh.write("(Sem NT — CDS real não mapeado.)\n")

    # regions.csv (resumo simples)
    rcsv = CONFIG["output"]["export_regions_csv"]
    if rcsv:
        try:
            os.makedirs(os.path.dirname(rcsv), exist_ok=True)
            with io.open(rcsv, "w", encoding="utf-8", newline="") as fcsv:
                w = csv.writer(fcsv)
                w.writerow(["Tag","Start","End","Length"])
                for a in aa_regions:
                    w.writerow([a["tag"], a["start"], a["end"], a["end"]-a["start"]+1])
        except Exception as e:
            print(f"[WARN] Falha ao escrever regions.csv: {e}")

    # variants_blast.csv (opcional; placeholder)
    vcsv = CONFIG["output"]["export_csv"]
    if vcsv and os.path.basename(vcsv):
        try:
            os.makedirs(os.path.dirname(vcsv), exist_ok=True)
            with io.open(vcsv, "w", encoding="utf-8", newline="") as fcsv:
                w = csv.writer(fcsv)
                w.writerow(["level","region_tag","pos","observed_set","source"])
        except Exception as e:
            print(f"[WARN] Falha ao escrever variants_blast.csv: {e}")

    # BLAST hits audit CSV (agora com species/same_species)
    hits_csv = CONFIG["output"].get("blast_hits_csv")
    if hits_csv:
        try:
            os.makedirs(os.path.dirname(hits_csv), exist_ok=True)
            with io.open(hits_csv, "w", encoding="utf-8", newline="") as fcsv:
                w = csv.writer(fcsv)
                w.writerow([
                    "kind","db","region","accession",
                    "percent_identity","query_coverage","hsp_query_len",
                    "species","same_species","definition"
                ])
                for rec in (_BLASTP_AUDIT + _BLASTN_AUDIT):
                    w.writerow([
                        rec.get("kind"), rec.get("db"), rec.get("region"), rec.get("accession"),
                        rec.get("percent_identity"), rec.get("query_coverage"), rec.get("hsp_query_len"),
                        rec.get("species",""), rec.get("same_species",""), rec.get("definition","")
                    ])
            if CONFIG["output"].get("blast_progress"):
                print(f"[OK] BLAST hits salvos em {hits_csv}")
        except Exception as e:
            print(f"[WARN] Falha ao escrever blast_hits_csv: {e}")

    # context_summary.json
    ctx_path = os.path.join(CONFIG["output"]["artifacts_dir"], "context_summary.json")
    try:
        os.makedirs(CONFIG["output"]["artifacts_dir"], exist_ok=True)
        with io.open(ctx_path, "w", encoding="utf-8") as fh:
            json.dump(context, fh, ensure_ascii=False, indent=2)
        print(f"[OK] context_summary.json salvo em {ctx_path}")
    except Exception as e:
        print(f"[ERRO] Ao salvar context_summary.json: {e}")

# =========================
# Entrypoint
# =========================
if __name__ == "__main__":
    # Sinalização amigável para NCBI email
    if not NCBI_EMAIL:
        print("[AVISO] Defina NCBI_EMAIL para cumprir a política do NCBI.")

    src = (CONFIG["input"]["source"] or "").lower()
    try:
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
    except requests.HTTPError as e:
        print(f"[ERRO HTTP] {e} — URL: {getattr(e.response, 'url', '?')}")
    except Exception as e:
        print(f"[ERRO] {e}")

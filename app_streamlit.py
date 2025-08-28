# app_streamlit.py
# MVP UI (Streamlit) — execução via SUBPROCESS (shim)
# - cria shim Python temporário que carrega seu main.py, aplica overrides no CONFIG e chama run_from_*.
# - captura stdout/stderr completos e mostra em "Logs".
# - lê artifacts/context_summary.json para renderizar cards/seqs/downloads.
# - mantém estado entre reruns; botão "Limpar" só zera entradas e resultados.
# - blocos de sequência usam st.code() (um só) com wrap visual via CSS e cópia sem quebras.
# - Gera automaticamente um arquivo "Prompt Effatha" (.md) com base no contexto + objetivo informado.

import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime

import streamlit as st

# -------------------------------
# CSS — wrap visual em blocos de código (sem inserir \n)
# -------------------------------
st.markdown(
    """
    <style>
    /* Força wrap visual em blocos de código (st.code),
       mas o conteúdo não contém quebras — copia sai contínuo */
    .stCode > div > pre, .stCode pre, pre, code {
        white-space: pre-wrap !important;
        word-break: break-word !important;
        overflow-wrap: anywhere !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# -------------------------------
# Utils
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
    for k, v in u.items():
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

# estado p/ manter último run
def save_session(run_dir, context, paths, enable_blast, prompt_goal, mech_flags):
    st.session_state["last_run_dir"] = str(run_dir)
    st.session_state["last_context"] = context
    st.session_state["last_paths"] = paths
    st.session_state["last_enable_blast"] = bool(enable_blast)
    # também lembrar últimas preferências de prompt
    st.session_state["last_prompt_goal"] = prompt_goal
    st.session_state["last_mech_flags"] = mech_flags

def clear_results_and_inputs():
    # Limpa apenas as entradas de fonte e os resultados
    for k in ("uniprot_acc","pdb_id","ncbi_protein_acc","gene_id","gene_symbol"):
        if k in st.session_state:
            st.session_state[k] = ""
    # Mantém gene_taxid e demais configs
    for k in ("last_run_dir","last_context","last_paths","last_enable_blast"):
        if k in st.session_state:
            st.session_state[k] = None

# -------------------------------------
# Prompt Effatha — geração de arquivo .md
# -------------------------------------
def _list_regions_lines(context):
    lines = []
    for r in context.get("regions", []) or []:
        note = f" — {r.get('note','')}" if r.get('note') else ""
        lines.append(f"- {r.get('source')}:{r.get('type')} {r.get('start')}-{r.get('end')}{note}")
    return "\\n".join(lines) if lines else "(sem regiões)"

def _seq_block(label, seq_text):
    if not seq_text:
        return ""
    return f"**{label}**\\n\\n```\n{seq_text}\n```\\n"

def _render_region_sequences(context):
    # prioriza OBS/FILTRADO; se ausentes, usa aa_regions (cru)
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
        blocks.append("\\n".join([x for x in section if x]))
    return "\\n\\n".join(blocks) if blocks else "_Sem sequências por região disponíveis._"

def build_effatha_prompt_md(*, context, treatment_goal, mech_flags, report_path):
    """
    Constrói o conteúdo Markdown do Prompt Effatha para uso em GPT Pro.
    mech_flags: dict com chaves 'aproximar', 'distanciar', 'mimetizar' (bool)
    """
    protein = context.get("protein", {}) or {}
    blast = context.get("blast", {}) or {}
    layers = context.get("layers", {}) or {}

    allow_list = [k for k,v in mech_flags.items() if v]
    allow_txt = ", ".join(allow_list) if allow_list else "aproximar, distanciar e mimetizar"

    region_lines = _list_regions_lines(context)
    region_seq_md = _render_region_sequences(context)

    report_tail = ""
    try:
        if report_path and Path(report_path).exists():
            txt = Path(report_path).read_text(encoding="utf-8", errors="ignore")
            # pega início do relatório para contexto (limita tamanho)
            head = txt.strip()
            if len(head) > 12000:
                head = head[:12000] + "\\n... (trecho truncado)"
            report_tail = f"\\n\\n---\\n### Trecho do relatório (.txt)\\n```\n{head}\n```"
    except Exception:
        pass

    md = f"""# Prompt Effatha — Análise e Priorização de Sequências (para GPT Pro)

**Objetivo do tratamento:** {treatment_goal or "(não informado)"}  
**Mecanismos permitidos nesta análise:** {allow_txt}

## 0) Papel da IA
Você é uma IA especialista em biologia molecular/sistemas, treinada para avaliar **potencial de intervenção por frequências Effatha** em proteínas/genes.  
Para **cada sequência** fornecida a seguir, você deve:
1. Classificar o **potencial** para Effatha como **"Alto"**, **"Moderado"**, **"Baixo"** ou **"Sem potencial"** (obrigatório considerar "Sem potencial").  
2. Indicar **qual mecanismo Effatha** é mais promissor (**aproximar**, **distanciar**, **mimetizar**) e por quê. Use apenas mecanismos **permitidos** acima.  
3. Descrever o **raciocínio** (máx. ~8 linhas por sequência), citando: função/região, contexto estrutural, variantes presentes (`[ref/alt/...]`), efeitos esperados, riscos/limitações.  
4. Retornar uma **priorização global** (Top N) com justificativas curtas, e **descartar** explicitamente alvos com “Sem potencial” (com uma frase de justificativa para economizar recursos).
5. Quando fizer sentido, sugerir **controles e próximos passos**.

> Observação: esta IA **não precisa** reproduzir exatamente o processo interno do Analyzer. Ela apenas precisa **entender** como as sequências foram derivadas (flancos, colchetes, etc.) para interpretar corretamente os blocos e realizar uma avaliação mais ampla.

## 1) Contexto do caso
- **Acesso**: {protein.get('accession','?')}
- **Organismo**: {protein.get('organism','?')} (taxid={protein.get('taxid','?')})
- **Tamanho (aa)**: {protein.get('length','?')}
- **Nome / Descrição**: {protein.get('name','')}
- **Camadas ativas (no Analyzer, para gerar contexto)**:
  - UniProt VARIANT = { "ON" if layers.get("uniprot_variant_enabled") else "OFF" }
  - Proteins Variation API = { "ON" if layers.get("proteins_variation_enabled") else "OFF" }
  - BLAST = { "ON" if layers.get("blast_enabled") else "OFF" }
- **BLAST (resumo)**: janelas={blast.get('windows',0)}, filtro={blast.get('filter','None')}, pos_var(bast)= {blast.get('variant_positions_blast',0)}

## 2) Como as sequências foram derivadas (Gene & Protein Analyzer)
- Regiões vêm de **UniProt/Pfam/NCBI** e podem incluir notas (ligantes, íons, domínios, **SITES** etc.).
- Cada região é expandida por **flancos** definidos (ex.: ±5 aa), e janelas são geradas/centralizadas.
- Sequências de AA podem conter **colchetes** com variantes ou alternativas, p.ex.: `A[L/E]K` (ref/alt/...).
- Também disponibilizamos **DNA** e **mRNA** correspondentes à região (de CDS real quando mapeado; senão, retrotradução).
- O objetivo **não** é repetir o pré-processamento, mas **interpretar** corretamente estes blocos.

## 3) Regiões-alvo (após expansão)
{region_lines}

## 4) Sequências por região (AA + DNA/mRNA)
{region_seq_md}

## 5) Mecanismos Effatha — definição operacional (resumo)
- **Aproximar**: reduzir distâncias efetivas/funcionais entre átomos/resíduos/regiões para **alterar conformação** (ex.: bloquear sítio ativo/ligação), **expor** ou **ocultar** elementos críticos, **inibir expressão gênica** (casos similares a RNAi observados em *T. rubrum*), ou **modular** interação PPI.  
- **Distanciar**: aumentar distâncias/afrouxar contatos para **expor sítios catalíticos**, **relaxar fibras** ou **facilitar interações** desejadas.  
- **Mimetizar**: simular frequência de ligantes/peptídeos/íon metálico para **induzir resposta semelhante** à presença química real **sem** adicionar a molécula.

## 6) Tarefa
Para **cada bloco de região** acima, retorne um item com:
- **Potencial Effatha**: Alto / Moderado / Baixo / **Sem potencial**.
- **Mecanismo(s) sugerido(s)** (restrito aos permitidos): Aproximar / Distanciar / Mimetizar.
- **Justificativa** (máx. ~8 linhas; cite função/ligante/domínio, variantes entre colchetes, efeitos esperados e riscos).
- **Teste/controle recomendado** (1–2 linhas).
Ao final, liste um **ranking dos melhores alvos** e uma **lista de descarte** (sem potencial) com 1 linha de motivo por item.

---
**Observações finais do operador**: Se necessário, considere a **não viabilidade experimental** (ex.: ausência de papel funcional claro, região fora do contexto da doença, variantes incompatíveis com o objetivo, etc.).{report_tail}
"""
    return md

def write_prompt_file(run_dir, *, context, treatment_goal, mech_flags, report_path, slug, ts):
    md = build_effatha_prompt_md(context=context,
                                 treatment_goal=treatment_goal,
                                 mech_flags=mech_flags,
                                 report_path=report_path)
    prompt_md_path = Path(run_dir) / f"prompt_effatha_{slug}_{ts}.md"
    prompt_txt_path = Path(run_dir) / f"prompt_effatha_{slug}_{ts}.txt"
    prompt_md_path.write_text(md, encoding="utf-8")
    prompt_txt_path.write_text(md, encoding="utf-8")
    return str(prompt_md_path), str(prompt_txt_path)

# -------------------------------
# UI — Sidebar / Parâmetros (com keys p/ podermos limpar)
# -------------------------------
st.set_page_config(page_title="Effatha — Functional Regions", layout="wide")
st.title("Gene & Protein Analyzer")

with st.sidebar:
    st.header("⚙️ Configurações")

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
        if st.button("Limpar"):
            clear_results_and_inputs()
            st.rerun()

    # Inputs por fonte (com keys para serem limpos)
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
    st.subheader("Saída")
    st.number_input("Largura de linha (visualização)", min_value=40, max_value=200, step=2, value=80, key="wrap_width")
    st.text_input("Pasta de artefatos", value="runs", key="artifacts_dir")

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

    btn = st.button("🚀 Executar análise", type="primary", key="run_btn")

tab_resumo, tab_regioes, tab_seqs, tab_logs, tab_downloads = st.tabs(
    ["Resumo", "Regiões", "Sequências", "Logs", "Downloads"]
)
log_placeholder = tab_logs.empty()
summary_placeholder = tab_resumo.empty()
regions_placeholder = tab_regioes.empty()
seqs_container = tab_seqs.container()
downloads_container = tab_downloads.container()

# -------------------------------
# Renderização
# -------------------------------
def render_outputs(context, paths, enable_blast_flag):
    # ----- Resumo -----
    protein = context.get("protein", {})
    layers = context.get("layers", {})
    blast = context.get("blast", {})
    with summary_placeholder:
        st.subheader("Resumo da Execução")
        colA, colB, colC = st.columns(3)
        with colA:
            st.metric("Acesso", protein.get("accession","?"))
            st.metric("Tamanho (aa)", protein.get("length","?"))
        with colB:
            st.metric("Espécie", f"{protein.get('organism','?')} (taxid={protein.get('taxid','?')})")
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

    # ----- Regiões (cards) -----
    cards = context.get("region_cards", []) or []
    if not cards:
        regions_placeholder.info("Sem region cards disponíveis.")
    else:
        import pandas as pd
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
            })
        df = pd.DataFrame(rows)
        regions_placeholder.dataframe(df, use_container_width=True)

    # ----- Sequências por região -----
    seqs_container.subheader("Sequências por região (AA com colchetes + DNA/mRNA)")
    aa_idx_obs = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions_obs", [])}
    aa_idx_fil = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions_filtered", [])}
    aa_idx_any = {(a["tag"], a["start"], a["end"]): a["aa"] for a in context.get("aa_regions", [])}
    cards = context.get("region_cards", []) or []
    for seg in context.get("nt_segments", []):
        tag = seg["tag"]; start = seg["start"]; end = seg["end"]
        c = pick_card(cards, tag, start, end)

        with st.expander(f"{tag} {start}-{end}"):
            # Cabeçalho + descrição
            header = f"**{tag}** {start}-{end}"
            if c and c.get("length"):
                header += f" (len={c.get('length')})"
            if seg.get("note"):
                header += f" — {seg['note']}"
            st.markdown(header)

            # Metadados
            if c:
                ligand_line = ""
                if c.get("ligand_name") or c.get("ligand_chebi"):
                    chebi_txt = f" ({c.get('ligand_chebi')})" if c.get("ligand_chebi") else ""
                    ligand_line = f"Ligante: {c.get('ligand_name')}{chebi_txt}"
                metal_line = ""
                if c.get("ligand_is_metal") is not None:
                    metal_line = f" · Íon metálico: {'sim' if c['ligand_is_metal'] else 'não'}"
                if ligand_line or metal_line:
                    st.caption((ligand_line or "") + (metal_line or ""))

                st.caption(
                    f"Variação: obs {c.get('pos_var_obs')} ({100.0*(c.get('pct_var_obs') or 0):.1f}%) · "
                    f"filtrado {c.get('pos_var_fil')} ({100.0*(c.get('pct_var_fil') or 0):.1f}%)"
                    + (" · deleção presente" if c.get("has_deletion") else "")
                )
                if c.get("max_af") is not None:
                    st.caption(f"AF máx (observado): {c.get('max_af'):.4f}")
                clin = c.get("clin_counts") or {}
                if clin:
                    st.caption("Clínica (observado): " + "; ".join([f"{k}={v}" for k,v in clin.items()]))
                st.caption(
                    f"Predições (observado): SIFT tol={c.get('sift_tol')}, del={c.get('sift_del')} | "
                    f"PolyPhen benign={c.get('poly_benign')}, damaging={c.get('poly_damaging')}"
                )

            # Sequências — APENAS um bloco (st.code) com wrap visual e cópia sem quebras
            aa_obs = aa_idx_obs.get((tag, start, end), "")
            aa_fil = aa_idx_fil.get((tag, start, end), "")
            aa_any = "" if (aa_obs or aa_fil) else aa_idx_any.get((tag, start, end), "")

            if aa_obs:
                st.markdown("**AA (observado)**")
                st.code(aa_obs)  # sem '\n' inseridos; wrap visual por CSS
            if aa_fil:
                st.markdown("**AA (filtrado)**")
                st.code(aa_fil)
            if aa_any:
                st.markdown("**AA**")
                st.code(aa_any)

            dna = seg.get("dna","")
            mrna = seg.get("mrna","")
            st.markdown("**DNA**")
            st.code(dna)
            st.markdown("**mRNA**")
            st.code(mrna)

    # ----- Downloads -----
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

            # NOVO: Prompt Effatha
            pmd = Path(paths.get("prompt_md","")) if paths else None
            ptxt = Path(paths.get("prompt_txt","")) if paths else None
            if pmd and pmd.exists():
                st.download_button("Baixar Prompt Effatha (.md)", data=pmd.read_bytes(),
                                   file_name=pmd.name, mime="text/markdown")
            if ptxt and ptxt.exists():
                st.download_button("Baixar Prompt Effatha (.txt)", data=ptxt.read_bytes(),
                                   file_name=ptxt.name, mime="text/plain")

        except Exception as e:
            st.error(f"Erro ao preparar downloads: {e}")

# -------------------------------
# Execução
# -------------------------------
if st.session_state.get("run_btn"):
    pp = Path(st.session_state["pipeline_path"])
    if not pp.is_file():
        st.error(f"Arquivo não encontrado: {pp}")
        st.stop()

    run_dir = make_run_dir(st.session_state.get("artifacts_dir") or "runs")

    # === nomes de arquivos com ID + timestamp ===
    def _sanitize(s: str) -> str:
        # mantém apenas letras/números e - _ .
        return "".join(ch for ch in (s or "") if ch.isalnum() or ch in ("-", "_", ".")).strip() or "run"

    src_label = st.session_state["source_choice"]
    if src_label == "UniProt":
        id_hint = st.session_state.get("uniprot_acc", "")
    elif src_label == "PDB":
        id_hint = st.session_state.get("pdb_id", "")
    elif src_label == "NCBI Protein":
        id_hint = st.session_state.get("ncbi_protein_acc", "")
    else:  # NCBI Gene
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
        "output": {
            "report_txt": report_path,
            "export_regions_csv": regions_csv_path,
            "export_csv": (blast_csv_path if st.session_state.get("enable_blast", False) else None),
            "artifacts_dir": str(run_dir),
        },
    }

    shim_path = write_shim(run_dir, st.session_state["pipeline_path"], cfg_overrides)
    py_exe = sys.executable

    with st.status("Executando pipeline…", expanded=True) as status:
        rc, out, err = run_subprocess(py_exe, shim_path, run_dir)
        log_placeholder.code((out or "") + ("\n" + err if err else ""))
        if rc == 0:
            status.update(label="Execução concluída", state="complete")
        else:
            status.update(label=f"Erro (exit code {rc}) — veja 'Logs'", state="error")

    context_path = Path(run_dir) / "context_summary.json"
    context = None
    if context_path.exists():
        try:
            context = json.loads(context_path.read_text(encoding="utf-8"))
        except Exception as e:
            tab_resumo.error(f"Falha ao ler context_summary.json: {e}")

    # NOVO: gerar Prompt Effatha automaticamente (se tivermos contexto)
    treatment_goal = st.session_state.get("treatment_goal","").strip()
    mech_flags = {
        "aproximar": bool(st.session_state.get("allow_aproximar", True)),
        "distanciar": bool(st.session_state.get("allow_distanciar", True)),
        "mimetizar": bool(st.session_state.get("allow_mimetizar", True)),
    }
    prompt_md_path = None
    prompt_txt_path = None
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
        render_outputs(context, paths, st.session_state.get("enable_blast", False))
    else:
        tab_resumo.info("Contexto não encontrado — verifique os Logs para erros no pipeline.")

elif st.session_state.get("last_context"):
    render_outputs(
        st.session_state.get("last_context"),
        st.session_state.get("last_paths", {}),
        st.session_state.get("last_enable_blast", False),
    )
else:
    summary_placeholder.info("Pronto para executar. Configure os parâmetros e clique em **Executar análise**.")

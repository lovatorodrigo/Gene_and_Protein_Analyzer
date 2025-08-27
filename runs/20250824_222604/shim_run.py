# AUTO-GERADO pelo app Streamlit (não editar)
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

module_path = r"""C:\Users\rodri\PycharmProjects\analyze_uniprot.py\main.py"""
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

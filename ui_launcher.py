# ui_launcher.py — v5 (1 aba garantida: headless forçado + detecção de servidor)
import os
import sys
import time
import socket
import subprocess
from pathlib import Path


def find_app_path(exe_dir: Path) -> Path:
    """
    Procura por app_streamlit.py em cenários comuns:
    - onedir: <exe_dir>/_internal/app_streamlit.py  (quando --add-data ...;_internal)
    - onedir (fallback incomum): <exe_dir>/_internal/_internal/app_streamlit.py
    - ao lado do exe: <exe_dir>/app_streamlit.py
    - onefile: <_MEIPASS>/app_streamlit.py         (quando --add-data ...;.)
    - onefile: <_MEIPASS>/_internal/app_streamlit.py (quando --add-data ...;_internal)
    """
    candidates = [
        exe_dir / "_internal" / "app_streamlit.py",
        exe_dir / "_internal" / "_internal" / "app_streamlit.py",
        exe_dir / "app_streamlit.py",
    ]
    for p in candidates:
        if p.exists():
            return p

    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        mp = Path(meipass)
        for p in [
            mp / "app_streamlit.py",
            mp / "_internal" / "app_streamlit.py",
        ]:
            if p.exists():
                return p

    raise FileNotFoundError(
        "Não encontrei 'app_streamlit.py'. Inclua com --add-data e garanta que caia em "
        "'_internal' (onedir) ou no root/_internal do _MEIPASS (onefile)."
    )


def write_streamlit_toml(exe_dir: Path, port: int = 8501) -> Path:
    """
    Gera um config TOML para o Streamlit:
    - headless = true  => Streamlit NÃO abre navegador sozinho
    - browser.serverAddress = ""  => inibe auto-open do Streamlit
    """
    cfg = (
        "[server]\n"
        "headless = true\n"
        f"port = {port}\n"
        "address = \"127.0.0.1\"\n"
        "fileWatcherType = \"none\"\n"
        "\n"
        "[global]\n"
        "developmentMode = false\n"
        "\n"
        "[browser]\n"
        "gatherUsageStats = false\n"
        "serverAddress = \"\"\n"
    )
    cfg_path = exe_dir / "effatha_streamlit.toml"
    cfg_path.write_text(cfg, encoding="utf-8")
    return cfg_path


def wait_for_port(port: int, host: str = "127.0.0.1", timeout_s: float = 15.0) -> bool:
    """Espera até a porta estar escutando (ou estoura timeout)."""
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(0.5)
            try:
                s.connect((host, port))
                return True
            except Exception:
                time.sleep(0.25)
    return False


def main() -> int:
    exe_dir = Path(sys.executable).resolve().parent
    log_path = exe_dir / "effatha_ui.log"
    port = 8501

    try:
        app_path = find_app_path(exe_dir)
        cfg_path = write_streamlit_toml(exe_dir, port=port)

        with open(log_path, "w", encoding="utf-8") as log:
            log.write(f"[INFO] Iniciando Streamlit: {app_path} (cfg={cfg_path})\n")

        # Força nosso TOML e headless (redundância para "calar" qualquer auto-open)
        os.environ["STREAMLIT_CONFIG"] = str(cfg_path)
        os.environ["STREAMLIT_SERVER_HEADLESS"] = "true"
        os.environ["STREAMLIT_GLOBAL_DEVELOPMENTMODE"] = "false"
        os.environ["STREAMLIT_BROWSER_SERVER_ADDRESS"] = ""
        # Desabilita qualquer tentativa de abrir navegador via módulo webbrowser
        os.environ["BROWSER"] = ""

        # Dispara o Streamlit via CLI programática
        # Reforçamos headless também por flag (para vencer qualquer config residual)
        py_snippet = (
            "import sys; import streamlit.web.cli as stcli; "
            "sys.argv = ['streamlit','run', r'%s', '--server.headless=true']; "
            "raise SystemExit(stcli.main())"
        ) % (str(app_path),)

        subprocess.Popen([sys.executable, "-c", py_snippet], cwd=str(exe_dir))

        # Espera o servidor subir antes de abrir a aba
        if wait_for_port(port, "127.0.0.1", timeout_s=20.0):
            try:
                import webbrowser
                webbrowser.open(f"http://localhost:{port}", new=2)
            except Exception:
                pass
        else:
            # Se não subiu no tempo esperado, registra no log
            with open(log_path, "a", encoding="utf-8") as log:
                log.write(f"[WARN] Servidor não respondeu na porta {port} dentro do timeout.\n")

        return 0

    except Exception as e:
        # Loga erro
        try:
            with open(log_path, "a", encoding="utf-8") as log:
                log.write(f"[ERRO] {e}\n")
                log.write("[TRACEBACK]\n")
                import traceback; log.write(traceback.format_exc())
        except Exception:
            pass

        # Caixa de diálogo no Windows (ignora em outros SOs)
        try:
            import ctypes  # type: ignore
            ctypes.windll.user32.MessageBoxW(
                0,
                f"Erro ao iniciar a UI:\n{e}\n\nVeja effatha_ui.log ao lado do executável.",
                "Effatha Protein Analyzer",
                0,
            )
        except Exception:
            pass

        return 1


if __name__ == "__main__":
    sys.exit(main())

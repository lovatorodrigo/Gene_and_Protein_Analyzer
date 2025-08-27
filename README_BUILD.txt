INSTRUÇÕES RÁPIDAS (ATUALIZADAS)

1) Copie estes arquivos para a RAIZ do seu projeto (junto de app_streamlit.py e main.py):
   - build_installer.bat
   - ui_launcher.py
   - effatha.iss
   - requirements.txt

2) Abra o Terminal do PyCharm (ou PowerShell) nessa pasta e rode:
   .\build_installer.bat

3) Saídas esperadas:
   - Portátil (onedir): dist\EffathaUI\EffathaUI.exe
   - Opcional onefile:  dist\Effatha\Effatha.exe
   - Instalador (se Inno Setup estiver instalado): Effatha-Setup.exe

Se ao abrir pelo atalho nada acontecer, veja o arquivo de log:
   dist\EffathaUI\effatha_ui.log

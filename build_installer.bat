@echo off
setlocal enableextensions enabledelayedexpansion

REM ===== Effatha Builder (limpo/estável) =====
REM Execute este .BAT na RAIZ do projeto (onde estao app_streamlit.py, main.py e ui_launcher.py).

cd /d "%~dp0"

echo.
echo [0/6] Checando Python...
where python >nul 2>nul
if errorlevel 1 (
  echo [ERRO] Python nao encontrado no PATH. Instale Python 3.9+ e tente novamente.
  pause
  exit /b 1
)

echo.
echo [1/6] Verificando .venv...
if not exist ".venv\Scripts\python.exe" (
  echo   Criando .venv...
  python -m venv .venv
  if errorlevel 1 (
    echo [ERRO] Falha ao criar ambiente virtual.
    pause
    exit /b 1
  )
) else (
  echo   .venv encontrada.
)

echo.
echo [2/6] Instalando dependencias...
call ".venv\Scripts\activate"
python -m pip install --upgrade pip
if exist "requirements.txt" (
  pip install -r requirements.txt
  if errorlevel 1 (
    echo [AVISO] Falha ao instalar requirements.txt. Instalando pacote basico.
    pip install streamlit requests biopython reportlab
  )
) else (
  echo [AVISO] requirements.txt nao encontrado. Instalando pacote basico.
  pip install streamlit requests biopython reportlab
)
pip install pyinstaller

echo.
echo [3/6] Limpando saidas anteriores...
if exist dist rd /s /q dist
if exist build rd /s /q build

echo.
echo [4/6] Gerando executavel ONEDIR...
set "APP_DIR=%CD%"
pyinstaller --noconfirm --name EffathaUI --onedir --noconsole ^
  --add-data "%APP_DIR%\app_streamlit.py;." ^
  --add-data "%APP_DIR%\main.py;." ^
  --collect-all streamlit ^
  --collect-all Bio ^
  ui_launcher.py

if errorlevel 1 (
  echo [ERRO] PyInstaller ONEDIR falhou. Verifique as mensagens acima.
) else (
  echo [OK] Portatil (onedir): dist\EffathaUI\EffathaUI.exe
)

echo.
echo [5/6] Gerando executavel ONEFILE (opcional)...
pyinstaller --noconfirm --name Effatha --onefile --noconsole ^
  --add-data "%APP_DIR%\app_streamlit.py;." ^
  --add-data "%APP_DIR%\main.py;." ^
  --collect-all streamlit ^
  --collect-all Bio ^
  ui_launcher.py

if errorlevel 1 (
  echo [AVISO] PyInstaller ONEFILE falhou. Ignorando e seguindo com onedir.
) else (
  echo [OK] Onefile: dist\Effatha\Effatha.exe
)

echo.
echo [6/6] Compilando instalador (Inno Setup), se disponivel...
REM --- descoberta do ISCC ---
set "ISCC="
if exist "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" set "ISCC=C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
if not defined ISCC if exist "C:\Program Files\Inno Setup 6\ISCC.exe" set "ISCC=C:\Program Files\Inno Setup 6\ISCC.exe"
if not defined ISCC (
  where ISCC.exe >nul 2>nul
  if not errorlevel 1 (
    for %%G in (ISCC.exe) do set "ISCC=%%~$PATH:G"
  )
)

if not defined ISCC goto NOINNO
if not exist "effatha.iss" goto NOISS

REM --- compila .iss usando caminhos do proprio .iss (usa SourcePath) ---
"%ISCC%" "effatha.iss"
if errorlevel 1 (
  echo [ATENCAO] Falha ao compilar o instalador via Inno Setup.
  goto THEEND
) else (
  echo [OK] Instalador gerado (procure em dist\installer\Effatha-Setup.exe).
  goto THEEND
)

:NOINNO
echo [AVISO] Inno Setup (ISCC.exe) nao encontrado. Instale para gerar Effatha-Setup.exe.
goto THEEND

:NOISS
echo [AVISO] effatha.iss nao encontrado na raiz do projeto. Pulando instalador.
goto THEEND

:THEEND
echo.
echo ===== FIM =====
echo - Portatil: dist\EffathaUI\EffathaUI.exe
echo - Onefile:  dist\Effatha\Effatha.exe
echo - Instalador: dist\installer\Effatha-Setup.exe (se Inno Setup estiver instalado)
pause

endlocal

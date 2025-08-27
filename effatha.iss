; effatha.iss — robusto com SourcePath
#define MyAppName        "Effatha Protein Analyzer"
#define MyAppVersion     "1.0"
#define AppPublisher     "Effatha"
#define AppExeName       "EffathaUI.exe"

#define AppRoot          SourcePath
#define AppBinDir        AppRoot + "\dist\EffathaUI"
#define OutputDirAbs     AppRoot + "\dist\installer"

[Setup]
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#AppPublisher}
DefaultDirName={pf}\{#MyAppName}
DefaultGroupName={#AppPublisher}
OutputDir={#OutputDirAbs}
OutputBaseFilename=Effatha-Setup
ArchitecturesInstallIn64BitMode=x64
DisableProgramGroupPage=yes
Compression=lzma2
SolidCompression=yes
WizardStyle=modern

[Files]
Source: "{#AppBinDir}\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#AppExeName}"
Name: "{commondesktop}\{#MyAppName}"; Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Criar atalho na área de trabalho"; GroupDescription: "Atalhos:"

[Run]
Filename: "{app}\{#AppExeName}"; Description: "Abrir agora"; Flags: nowait postinstall skipifsilent

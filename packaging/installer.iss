; Inno Setup script for the qPCR Analyzer Windows installer.
; Builds a standard "Next-Next-Finish" installer that drops the
; PyInstaller one-folder build under %LOCALAPPDATA%\Programs and adds
; Start-Menu + optional desktop shortcuts.
;
; Compiled in CI by Inno Setup 6 (ISCC.exe). The CI workflow passes the
; version via /DAppVersion=X.Y.Z so it stays in sync with pyproject.toml.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif

#define AppName        "qPCR Analyzer"
#define AppPublisher   "Jielin Yang"
#define AppURL         "https://github.com/j-y26/py_qpcr_analyzer"
#define AppExeName     "qpcr-analyzer.exe"
#define SourceDir      "..\dist\qpcr-analyzer"

[Setup]
AppId={{A6F2C0C2-3F2A-4E5C-9B8F-7B0F4A4E2D11}}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
AppSupportURL={#AppURL}/issues
AppUpdatesURL={#AppURL}/releases
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
LicenseFile=..\LICENSE
OutputDir=..\dist\installer
OutputBaseFilename=qpcr-analyzer-{#AppVersion}-windows-x64-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
UninstallDisplayIcon={app}\{#AppExeName}
UninstallDisplayName={#AppName} {#AppVersion}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExeName}"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExeName}"; Description: "Launch {#AppName}"; Flags: nowait postinstall skipifsilent

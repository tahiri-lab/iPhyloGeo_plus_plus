#define MyAppName "iPhyloGeo++"
#define MyAppPublisher "tahirilab"
#define MyAppURL "https://github.com/tahiri-lab/iPhyloGeo_plus_plus"
#define MyAppExeName "iPhyloGeo_plus_plus.exe"

[Setup]
AppName=iPhyloGeo++ and Ghostscript
AppVersion=1.0
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={pf}\tahirilab\{#MyAppName}
DefaultGroupName=iPhyloGeo++ and Ghostscript Group
LicenseFile=..\..\exe-generation\license.txt
Uninstallable=no
CreateAppDir=no
OutputBaseFilename=iPhyloGeo++andGhostscript_bundledinstaller

[Files]
Source: "gs10060w64.exe"; DestDir: "{tmp}"
Source: "iPhyloGeo++-1.0-win64.msi"; DestDir: "{tmp}"

[Run]
Filename: "{tmp}\gs10060w64.exe"; Description: "Install Ghostscript (Artifex Software)"; Flags: skipifsilent
Filename: "msiexec.exe"; Parameters: "/i ""{tmp}\iPhyloGeo++-1.0-win64.msi"""
[Setup]
AppName=iPhyloGeo++ and Ghostscript
AppVersion=1.0
DefaultDirName={pf}\tahirilab\iPhyloGeo_plus_plus_and_Ghostscript
DefaultGroupName=iPhyloGeo++ and Ghostscript Group
Uninstallable=no
CreateAppDir=no

[Files]
Source: "gs10060w64.exe"; DestDir: "{tmp}"
Source: "iPhyloGeo++-1.0-win64.msi"; DestDir: "{tmp}"

[Run]
Filename: "{tmp}\gs10060w64.exe"; Description: "Install Ghostscript (Artifex Software)"; Flags: skipifsilent
Filename: "msiexec.exe"; Parameters: "/i ""{tmp}\iPhyloGeo++-1.0-win64.msi"""
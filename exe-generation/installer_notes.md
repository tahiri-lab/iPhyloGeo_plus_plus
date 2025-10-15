# Installer Creation Notes

## October 13, 2025

1. Downloaded [Inno Setup](https://jrsoftware.org/isdl.php)
2. Put the app licenses in a txt file to have the MSI display them
3. Used the Inno Setup wizard to create an initial .iss
4. Downloaded the 64-bit Windows versions of Ghostscript
5. In the .iss, added `Source: "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\exe-generation\gs10060w64.exe"; DestDir: "{tmp}"; Flags: ignoreversion` to `[Files]`
6. In the .iss, at the very start of `[Run]`, added `Filename: "{tmp}\gs10060w64.exe"; Description: "Install Ghostscript (Artifex Software)"; Flags: postinstall skipifsilent`

## Octobre 15, 2025

1. Created an installer an rand it (chose to install only for the active user, and not to install Ghostscript since I already had it installed)
2. When trying to launch the freshly installed iPhyloGeo++, I get a **cx_Freeze fatal error**: `failed to get the Python codec of the filesystem encoding`
3. Opened my AppData\Local\Programs\iPhyloGeo++ folder, everything seemed OK
4. Deleted the folder and the shortcut from the start menu

HT create an installer
1. Install Inno Setup
2. Download the latest 64-bit Windows version of Ghostscript
https://www.ghostscript.com/releases/gsdnld.html
3. Open iPhyloGeo_plus_plus.iss in Inno Setup
4. In the `[Files]` section, replace the iPhyloGeo_plus_plus directory’s path with the correct path to it on your machine
5. In the `[Files]` section, replace the path to the Ghostscript installer (second line in `[Files]`) with the correct path to it on your machine
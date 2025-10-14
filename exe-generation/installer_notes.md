# Installer Creation Notes

## October 13, 2025

1. Downloaded [Inno Setup](https://jrsoftware.org/isdl.php)
2. Put the app licenses in a txt file to have the MSI display them
3. Created a first version of the installer
4. Downloaded both Windows versions of Ghostscript
5. In the .iss, added `Source: "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\exe-generation\gs10060w64.exe"; DestDir: "{tmp}"; Flags: ignoreversion` to `[Files]`
6. In the .iss, at the very start of `[Run]`, added `Filename: "{tmp}\gs10060w64.exe"; Description: "Install Ghostscript (Artifex Software)"; Flags: postinstall skipifsilent`


HT create an installer
1. Install Inno Setup
2. Download the latest 64-bit Windows version of Ghostscript
https://www.ghostscript.com/releases/gsdnld.html
3. Open iPhyloGeo_plus_plus.iss in Inno Setup
4. In the `[Files]` section, replace the iPhyloGeo_plus_plus directory’s path with the correct path to it on your machine
5. In the `[Files]` section, replace the path to the Ghostscript installer (second line in `[Files]`) with the correct path to it on your machine
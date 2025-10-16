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
5. In scripts/setup.py, changed base from "gui" to "Win32GUI".
6. Generated a new iPhyloGeo++ EXE
7. Generated a new installer EXE and ran it: upon running iPhyloGeo++.exe, I got a new error about the **DDL search path**: `Unable to change DLL search path`
8. Tried running the EXE as administrator: same error
8. Researched the error
9. Researched how Windows handles DLLs
10. Researched [packaged and unpackaged apps](https://learn.microsoft.com/en-us/windows/apps/get-started/intro-pack-dep-proc). Packaged apps seem to only be installable at a user level, not a system-wide level. Thus iPhyloGeo++ should remain unpackaged.
11. Further researched [how Windows handles DLLs](https://learn.microsoft.com/en-us/windows/win32/dlls/dynamic-link-library-search-order)
12. Researched [DLLs](https://learn.microsoft.com/en-ca/troubleshoot/windows-client/setup-upgrade-and-drivers/dynamic-link-library)
13. Created a [Stack Overflow question](https://stackoverflow.com/questions/79791613/cx-freeze-fatal-error-unable-to-change-dll-search-path) to ask for help
14. Checked the EXE’s folder: there is nothing there that isn’t included in the .iss `[Files]` section
15. Rebuilt the frozen app because I remember a warning from cx_Freeze at the end. Found the warning, which listed multiple possible missing dependencies.

Missing dependencies:
? LIBPQ.dll
? MIMAPI64.dll
? OCI.dll
? Qt63DQuickScene3D.dll
? Qt6VirtualKeyboardQml.dll
? WINSPOOL.DRV
? api-ms-win-core-heap-l2-1-0.dll
? api-ms-win-core-libraryloader-l1-2-0.dll
? api-ms-win-core-libraryloader-l1-2-1.dll
? api-ms-win-core-path-l1-1-0.dll
? api-ms-win-core-realtime-l1-1-1.dll
? api-ms-win-core-winrt-l1-1-0.dll
? api-ms-win-core-winrt-string-l1-1-0.dll
? api-ms-win-power-base-l1-1-0.dll
? api-ms-win-shcore-scaling-l1-1-1.dll
? bthprops.cpl
? fbclient.dll
This is not necessarily a problem - the dependencies may not be needed on this platform.


HT create an installer
1. Install Inno Setup
2. Download the latest 64-bit Windows version of Ghostscript
https://www.ghostscript.com/releases/gsdnld.html
3. Open iPhyloGeo_plus_plus.iss in Inno Setup
4. In the `[Files]` section, replace the iPhyloGeo_plus_plus directory’s path with the correct path to it on your machine
5. In the `[Files]` section, replace the path to the Ghostscript installer (second line in `[Files]`) with the correct path to it on your machine
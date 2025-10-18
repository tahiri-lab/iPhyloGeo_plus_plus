# Installer Creation Notes

## October 13, 2025

1. Downloaded [Inno Setup](https://jrsoftware.org/isdl.php)
2. Put the app licenses in a txt file to have the MSI display them
3. Used the Inno Setup wizard to create an initial .iss
4. Downloaded the 64-bit Windows versions of Ghostscript
5. In the .iss, added `Source: "C:\Users\agaco\Documents\Phylogeo\iPhyloGeo_plus_plus\exe-generation\gs10060w64.exe"; DestDir: "{tmp}"; Flags: ignoreversion` to `[Files]`
6. In the .iss, at the very start of `[Run]`, added `Filename: "{tmp}\gs10060w64.exe"; Description: "Install Ghostscript (Artifex Software)"; Flags: postinstall skipifsilent`

## October 15, 2025

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

## October 16, 2025

1. Searched my entire hard drive for the first possible missing dependency on the list, libpq.dll. It doesn’t exist in the project or its venv, but I do have several different versions. If needed, there is one in `C:\Users\myusername\anaconda3\pkgs\libpq-12.17-h906ac69_0\Library\bin\libpq.dll`

## October 17, 2025

New strategy: manually performing the steps I expect the installer to perform, on Windows Sandbox.

1. Created a ZIP archive of the exe.win-amd64-3.11 directory
2. Renamed the archive to iPhyloGeo_plus_plus
3. Copied both the iPhyloGeo_plus_plus.zip and the Ghostscript installer to Windows Sandbox
4. Unpacked the archive in C:\Program Files\iPhyloGeo_plus_plus in Windows Sandbox
5. Installed Ghostscript in Windows Sandbox
6. Ran iPhyloGeo_plus_plus.exe in Windows Sandbox: **it works!**
7. Reset Windows Sandbox
8. Copied my app’s installer to Windows Sandbox
9. Ran the installer, choosing a system-wide install and running the Ghostscript install
10. Ran the app: **cx_Freeze Fatal Error**: "Unable to change DLL search path"
11. Reset Windows Sandbox and copied the installer to it
11. Ran the installer, choosing a user install and running the Ghostscript install
12. Ran the app: **cx_Freeze Fatal Error**"Unable to change DLL search path"

New strategy: creating an MSI installer using cx_Freeze itself

1. Edited setup.py to add a bdist_msi_options with a randomly generated GUID as an upgrade code and set the initial_target_dir to `rf'[ProgramFilesFolder]\tahirilab\{app_name}'`
2. In the build .bat script, replaced `poetry run python setup.py build` with `poetry run python setup.py bdist_msi`
3. Ran the build .bat script
4. Launched Windows Sandbox
5. Copied the Ghostscript installer to Windows Sandbox
6. Ran the Ghostscript install
7. Looked at the results of the build: **ValueError**: `upgrade-code must be in valid GUID format`. The upgrade code was a 128-bit one, `'3216444d-fd6f-4b0a-a7e6-3b8b33ac5e5'`, so only HEX notation and dashes; it doesn’t have brackets.
8. Changed the upgrade code to `'{3216444d-fd6f-4b0a-a7e6-3b8b33ac5e54}'` (added curly brackets) and ran the build again
9. Copied the MSI to Windows Sandbox
10. Ran the MSI in Windows Sandbox
11. Launched the app: **success**

I still need to create an installer that will bundle the Ghostscript installer, and add a shortcut to the start menu.

Very useful: [https://ggottemo.com/blog/CxFreeze]

## October 18, 2025

1. Downloaded the .NET 9.0 SDK installer in order to get the WiX toolset from Microsoft.
2. Installed the .NET SDK

I decided to first follow a tutorial to familiarize myself with WiX: [https://docs.firegiant.com/quick-start/] Then, I got into actually creating a bundled installer for iPhyloGeo++ and Ghostscript

1. Created exe-generation\bundled_installer_files to store the necessary files to create a bundled Installer
2. Created exe-generation\bundled_installer_files\iPhyloGeo++.wixproj, a simple file indicating the SDK to use
3. Created exe-generation\bundled_installer_files\Package.wxs, with the contents of the first code snippet from [this Stack Overflow answer](https://stackoverflow.com/a/42102377/8814975)



HT create an installer that bundles Ghostscript
1. If you don’t already have it, get the [.NET SDK](https://dotnet.microsoft.com/en-us/download)
2. Build the project’s EXE and MSI using exe-generation\BUILDwMSI_frozenapp.bat

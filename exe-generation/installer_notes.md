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
4. Edited Package.wxs, setting the Bundle Name to "iPhyloGeo++ and Ghostscript", the UpgradeCode to the one from the setup.py, the LicenceUrl to [https://opensource.org/license/MIT], the SourceFile of the MSI to "iPhyloGeo++-1.0-win64.msi" and the SourceFile for the EXE to "gs10060w64.exe". Also made up package ids (text).
5. Moved the Ghostscript installer to the bundled_installer_files directory
6. Copied the contents of bundled_installer_files to scripts\dist
7. Opened a terminal from scripts\dist and ran `dotnet build`: got a **WIX0199 error**: Incorect namespace.
8. Swapped the url in line 1 of Package.wxs with "http://wixtoolset.org/schemas/v4/wxs"
9. Tried the build again: got two **errors**: **WIX0200** and **WIX0414**
10. Examined the details of error WIX0414: it says that the ExePackage element is missing a 'DetectCondition' attribute or a 'ArpEntry' child element because the 'Permanent' attribute is not specified
11. Added the Permanent attribute to the tag, set to "yes"
12. Examined the details of error WIX0200: "The BootstrapperApplicationRef element contains an unhandled extension element 'WixStandardBootstrapperApplication'. Please ensure that the extension for elements in the 'http://schemas.microsoft.com/wix/BalExtension' namespace has been provided". Decided to go for a more recent version found [here](https://docs.firegiant.com/wix/tools/burn/wixstdba/)
13. Tried the build again: WIX0200 still happens, whith a similar message as before. I also get a WIX0044 error: "The BootstrapperApplication element's Name or SourceFile attribute was not found; one of these is required."
14. Commented out the BootstrapperApplication element from Package.wxs
15. Tried the build again: **error WIX0350**: "The package being validated requires a higher version of Windows Installer than is installed on this machine. Validation cannot continue."
16. Researched the error and only found unanswered questions on GitHub and StackOverflow

## October 23, 2025

1. Read [https://docs.firegiant.com/wix/tools/validation/]
2. In iPhyloGeo++.wixproj, added a `PropertyGroup` tag containing a `SuppressValidation` tag with contents `true`
3. In the dist folder, ran `dotnet build`. It displayed 17 warnings but encountered no error
4. There was a new installer in bin\Debug, but it was extremely small (20 KB). Trying to run it triggered an **error message** saying I had to **install a Windows service pack containing a more recent version of Windows Installer**.
5. Used Window + R to run `MsiExec`: I have Windows ® Installer version 5.0.26100.5074, which is not the most recent version
6. Downloaded Windows 11, version 25H2
7. Installed Windows 11, version 25H2
8. Went back to scripts\dist, deleted the subfolders and ran `dotnet build`: same error message as before the update
9. Checked Windows Installer version: still 5.0.26100.5074
10. Checked Windows Update settings. The only optional updates it offers to install are drivers for my hardware
11. It doesn’t appear that there is a more recent version of Windows Installer I can install [https://learn.microsoft.com/en-ca/windows/win32/msi/windows-installer-redistributables]
12. Researched the error 'This installation package cannot be installed by the Windows Installer service. You must install a Windows service pack that contains a newer version of the Windows Installer service'
13. Ran in command line as administrator: `msiexec.exe /unregister` then `msiexec.exe /regserver`
14. Determined that the error message being incorrect was likely linked to me disabling errors in WiX. Commented out the lines added at step 2. I am now back to getting Error Wix0350 and, now that I’m reading it, the error message still makes no sense.
15. Restarted my machine
16. Downloaded the Ghostscript installer again
17. Commented out the EXE install from the Package.wxs file to make sure the issue is in fact with Ghostscript: same error, so it’s my installer (created with cx_Freeze) it doesn’t like
18. Reread the Package.wxs file and noted that I came up with the "ExePackage Id" values without following any rules as I was lacking information. Maybe there is a specific way to come up with IDs and my failure to follow it resulted in the issue?
19. Found a definition [here](https://docs.firegiant.com/wix3/xsd/wix/exepackage/): "Identifier for this package, for ordering and cross-referencing. The default is the Name attribute modified to be suitable as an identifier (i.e. invalid characters are replaced with underscores)." The id is a string. The documentation does not specify which characters are considered valid. The values I came up with only contained letters and underscores.
20. In Package.wxs, put the EXE install back in and commented out the MSI install instead: the error still occured
21. Drafted a Stackoverflow question in case Inno Setup fails to solve this

I decided to try Inno Setup instead.

1. Started from the code snippet from [this Stack Overflow answer](https://stackoverflow.com/a/15746747/) as a .iss
2. In the .iss, edited the files section to contain the actual names of the EXE and MSI
3. In the .iss, changed the first line in the `[Run]` section to `Filename: "{tmp}\gs10060w64.exe"; Description: "Install Ghostscript (Artifex Software)"; Flags: skipifsilent`
4. In the .iss, edited the last line to contain the name of the MSI
5. Ran the build
6. Launched Windows Sandbox and copied mysetup.exe to it
7. In Windows Sandbox, ran mysetup.exe: it unpacked both installers correctly, then ran the Ghostscript installer, then ran the iPhyloGeo++ installer, as expected

I later asked a friend who didn’t have Python or Ghostscript installed to test the bundled installer on his Windows 11 machine. He succeeded and was able to use iPhyloGeo++.

## October 25

1. Consulted [this StackOverflow answer](https://stackoverflow.com/a/15736406/) and [this table](https://learn.microsoft.com/en-us/windows/win32/msi/property-reference#system-folder-properties)
2. In scripts\setup.py, set shortcut_name to `app_name` and shortcut_dir to `"StartMenuFolder"` so that iPhyloGeo++ is added to the user’s start menu by the iPhyloGeo++ installer
3. Ran the build (with MSI)
4. In the .iss file’s Setup section, set the OutputBaseFilename variable in order to give a specific name to the generated installer
5. Built the bundled installer

TODO 

Check if I need to add a license (ask Nadia)
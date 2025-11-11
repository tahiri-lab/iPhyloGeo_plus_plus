from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.

dev_name = "tahirilab"

app_name = "iPhyloGeo++"

# Generate a MSI installer
bdist_msi_options = {
    'upgrade_code': '{3216444d-fd6f-4b0a-a7e6-3b8b33ac5e54}', # this is a randomly generated GUID; do not modify it so future versions of iPhyloGeo++ are recognized as the same app
    'add_to_path': False,
    'initial_target_dir': r'[ProgramFilesFolder]\%s\%s' % (dev_name, app_name),
    }

# Ensure icons from the img\disabled directory are included.
to_include = ["climatic", "climaticData", "download", "genetic", "result", "sequence", "setbutton", "settings", "start", "statistics", "submit", "tree"]
includefiles = []
for image_name in to_include:
    includefiles.append((f'..\img\disabled\{image_name}.svg', f'img\disabled\{image_name}.svg'))

build_options = {'packages': [], 'excludes': [], 'include_files':includefiles, 'optimize':2}

base = 'Win32GUI'

executables = [
    Executable('main.py', base=base, target_name= f'{app_name}.exe', shortcut_name=app_name, shortcut_dir="StartMenuFolder")
]

setup(name=app_name,
      version = '1.0',
      description = 'Multi-platform Application for Analyzing Phylogenetic Trees with Climatic Parameters',
      options = {'bdist_msi': bdist_msi_options, 'build_exe': build_options},
      executables = executables)

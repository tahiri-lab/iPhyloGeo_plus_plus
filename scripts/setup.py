from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.

app_name = "iPhyloGeo++"

# Generate a MSI installer
bdist_msi_options = {
    'upgrade_code': '3216444d-fd6f-4b0a-a7e6-3b8b33ac5e54', # this is a randomly generated GUID; do not modify it so future versions of iPhyloGeo++ are recognized as the same app
    'initial_target_dir': rf'[ProgramFilesFolder]\tahirilab\{app_name}'
    }

# Ensure icons from the img\disabled directory are included.
to_include = ["climatic", "climaticData", "download", "genetic", "result", "sequence", "setbutton", "settings", "start", "statistics", "submit", "tree"]
includefiles = []
for image_name in to_include:
    includefiles.append((f'..\img\disabled\{image_name}.svg', f'img\disabled\{image_name}.svg'))

build_options = {'packages': [], 'excludes': [], 'include_files':includefiles, 'optimize':1}

base = 'Win32GUI'

executables = [
    Executable('main.py', base=base, target_name= f'{app_name}.exe')
]

setup(name=app_name,
      version = '1.0',
      description = 'Multi-platform Application for Analyzing Phylogenetic Trees with Climatic Parameters',
      options = {'bdist_msi': bdist_msi_options, 'build_exe': build_options},
      executables = executables)

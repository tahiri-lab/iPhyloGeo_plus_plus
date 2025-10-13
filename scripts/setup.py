from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.

# Ensure icons from the img\disabled directory are included.
to_include = ["climatic", "climaticData", "download", "genetic", "result", "sequence", "setbutton", "settings", "start", "statistics", "submit", "tree"]
includefiles = []
for image_name in to_include:
    includefiles.append((f'..\img\disabled\{image_name}.svg', f'img\disabled\{image_name}.svg'))

build_options = {'packages': [], 'excludes': [], 'include_files':includefiles}

base = 'gui'

executables = [
    Executable('main.py', base=base, target_name='iPhyloGeo_plus_plus.exe')
]

setup(name='iPhyloGeo_plus_plus',
      version = '1.0',
      description = 'Multi-platform Application for Analyzing Phylogenetic Trees with Climatic Parameters',
      options = {'build_exe': build_options},
      executables = executables,
      optimize = 2)

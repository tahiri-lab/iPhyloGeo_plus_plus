from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
build_options = {'packages': [], 'excludes': []}

base = 'gui'

executables = [
    Executable('main.py', base=base)
]

setup(name='iPhyloGeo_plus_plus',
      version = '1.0',
      description = 'Multi-platform Application for Analyzing Phylogenetic Trees with Climatic Parameters',
      options = {'build_exe': build_options},
      executables = executables)

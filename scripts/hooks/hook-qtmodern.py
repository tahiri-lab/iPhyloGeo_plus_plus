from PyInstaller.utils.hooks import collect_data_files

datas = collect_data_files("qtmodern", includes=["**/*.qss"])

hiddenimports = [
    "qtmodern",
    "qtmodern.styles",
    "qtmodern.windows",
]
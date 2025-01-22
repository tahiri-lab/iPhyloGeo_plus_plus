if ($null -eq (Get-Command "pyuic5" -ErrorAction SilentlyContinue))
{
  Write-Host "Unable to find pyuic5 in PATH. Please install it by activating the venv using 'pip install pyqt5-tools'."
} else
{
  pyuic5 --import-from=utils -x .\scripts\Qt\main.ui -o .\scripts\Qt\main.py
  pyuic5 --import-from=utils -x .\scripts\Qt\loading.ui -o .\scripts\Qt\loading.py
  pyuic5 --import-from=utils -x .\scripts\Qt\parameters.ui -o .\scripts\Qt\parameters.py
  pyuic5 --import-from=utils -x .\scripts\Qt\wait.ui -o .\scripts\Qt\wait.py
  pyuic5 --import-from=utils -x .\scripts\Qt\help.ui -o .\scripts\Qt\help.py
}

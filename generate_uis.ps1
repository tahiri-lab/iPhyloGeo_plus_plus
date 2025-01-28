if ($null -eq (Get-Command "pyuic6" -ErrorAction SilentlyContinue))
{
  Write-Host "Unable to find pyuic6 in PATH. Please install it by activating the venv using 'pip install pyqt5-tools'."
} else
{
  pyuic6 -x .\scripts\Qt\main.ui -o .\scripts\Qt\main.py
  pyuic6 -x .\scripts\Qt\loading.ui -o .\scripts\Qt\loading.py
  pyuic6 -x .\scripts\Qt\parameters.ui -o .\scripts\Qt\parameters.py
  pyuic6 -x .\scripts\Qt\wait.ui -o .\scripts\Qt\wait.py
  pyuic6 -x .\scripts\Qt\help.ui -o .\scripts\Qt\help.py
}

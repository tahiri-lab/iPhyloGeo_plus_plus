if ($null -eq (Get-Command "pyuic6" -ErrorAction SilentlyContinue))
{
  Write-Host "Unable to find pyuic6 in PATH. Please install it by activating the venv using 'pip install pyqt5-tools'."
} else
{
  pyuic6 -x .\scripts\Qt\main.ui -o .\scripts\Qt\main_ui.py
  pyuic6 -x .\scripts\Qt\loading.ui -o .\scripts\Qt\loading_ui.py
  pyuic6 -x .\scripts\Qt\parameters.ui -o .\scripts\Qt\parameters_ui.py
  pyuic6 -x .\scripts\Qt\wait.ui -o .\scripts\Qt\wait_ui.py
  pyuic6 -x .\scripts\Qt\help.ui -o .\scripts\Qt\help_ui.py
}

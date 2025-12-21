import os
import subprocess
from pathlib import Path
import pytest

from aphylogeo import utils
from aphylogeo.alignement import Alignment
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
from aphylogeo.utils import fasttree


class TestFunctionalClimate:
    """
    Tests fonctionnels de l'exécutable iPhyloGeo++ pour la partie climatique.
    """

    def setup_class(self):
        self.base_dir = Path(__file__).resolve().parent.parent
        self.exe_path = self.base_dir / "scripts" / "build" / "exe.win-amd64-3.12" / "iPhyloGeo++.exe" 
        self.input_file = self.base_dir / "datasets" / "geo_with_loc.csv"
        assert self.exe_path.exists(), f"Executable not found: {self.exe_path}"

    def test_climate_analysis_execution(self):
        """
        Vérifier que l'exécutable s’exécute sans erreur et produit une sortie valide.
        """
        result = subprocess.run(
            [str(self.exe_path), "--mode", "climate", "--input", str(self.input_file)],
            capture_output=True, text=True
        )

        assert result.returncode == 0, f"Execution error: return code {result.returncode} {result.stderr}"

import os
import subprocess
import pytest
from pathlib import Path


from aphylogeo import utils
from aphylogeo.alignement import Alignment
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
from aphylogeo.utils import fasttree


class TestFunctionalGenetic:
    """
    Test de l'exécutable iPhyloGeo++ pour la partie génétique.
    """

    def setup_class(self):
        self.base_dir = Path(__file__).resolve().parent.parent
        self.exe_path = self.base_dir / "scripts" / "build" / "exe.win-amd64.3.12" / "iPhyloGeo++.exe" 
        self.input_file = self.base_dir / "datasets" / "example" / "geo.csv"
        assert self.exe_path.exists(), f"Executable non trouvé : {self.exe_path}"

    def test_genetic_analysis_execution(self):
        """
        Vérifier que l'exécutable s’exécute sans erreur et produit une sortie valide.
        """
        result = subprocess.run(
            [str(self.exe_path), "--mode", "genetic", "--input", str(self.input_file)],
            capture_output=True, text=True
        )

        # Vérifier que le processus se termine correctement
        assert result.returncode == 0, f"Erreur d'exécution : {result.stderr}"
        # Vérifier que la sortie contient une indication de succès
        assert "Analyse génétique terminée" in result.stdout

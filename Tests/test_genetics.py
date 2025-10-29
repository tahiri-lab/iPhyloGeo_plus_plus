import os
import subprocess
from pathlib import Path
import pytest

from aphylogeo import utils
from aphylogeo.alignement import Alignment
from aphylogeo.alignement import AlignSequences
from aphylogeo.params import Params
from aphylogeo.utils import fasttree

BASE_DIR = Path(__file__).resolve().parents[1]

EXECUTABLE = BASE_DIR / "exe-generation" / "iphylogeo++.exe"

TEST_DIR = BASE_DIR / "Test"

OUTPUT_DIR = TEST_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

class TestGeneticExecutable:

    def test_genetic_analysis_execution(self):

        print(">>> Test fonctionnel : exécution de l'analyse génétique")

        input_file = TEST_DIR / "datasets" / "gene_sequences.fasta"

        output_file = OUTPUT_DIR / "result_genetic.csv"

        cmd = [
            str(EXECUTABLE),
            "--mode", "genetic",
            "--input", str(input_file),
            "--output", str(output_file)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        assert result.returncode == 0, f"Erreur d’exécution : {result.stderr}"

        assert output_file.exists(), "Le fichier de sortie génétique n’a pas été généré."

        assert output_file.stat().st_size > 0, "Le fichier de sortie génétique est vide."


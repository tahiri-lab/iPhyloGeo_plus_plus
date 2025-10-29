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


DATASETS_DIR = TEST_DIR / "datasets" / "example"


OUTPUT_DIR = TEST_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def test_climatic_pipeline_functional():

    print(">>> Début du test fonctionnel : pipeline climatique")

    input_file = DATASETS_DIR / "geo.csv"

    output_file = OUTPUT_DIR / "climatic_result.json"

    cmd = [
        str(EXECUTABLE),
        "--mode", "climate",
        "--input", str(input_file),
        "--output", str(output_file)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"Erreur d'exécution : {result.stderr}"

    assert output_file.exists(), "Le fichier de sortie n'a pas été généré."

    assert output_file.stat().st_size > 0, "Le fichier de sortie est vide."

    print(">>> Test fonctionnel du pipeline climatique terminé avec succès.")


def test_climatic_dissimilarity_and_tree():

    print(">>> Début du test fonctionnel : dissimilarités et arbres climatiques")
    
    input_file = DATASETS_DIR / "geo.csv"

    output_file = OUTPUT_DIR / "climatic_tree.txt"

    cmd = [
        str(EXECUTABLE),
        "--mode", "climate-tree",
        "--input", str(input_file),
        "--output", str(output_file)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"Erreur d'exécution : {result.stderr}"
    assert output_file.exists(), "Le fichier d'arbre climatique n'a pas été généré."
    assert output_file.stat().st_size > 0, "Le fichier d'arbre climatique est vide."

    print(">>> Test fonctionnel des arbres climatiques terminé avec succès.")

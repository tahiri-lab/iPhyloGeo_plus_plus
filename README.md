<h1 align="center">🌳 iPhyloGeo++</h1>

<p align="center">
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
  <a href="https://pysd.readthedocs.io/en/latest/development/development_index.html">
    <img src="https://img.shields.io/badge/Contributions-Welcome-blue.svg" alt="Contributions">
  </a>
  <a href="https://pypi.python.org/pypi/pysd/">
    <img src="https://img.shields.io/pypi/pyversions/pysd.svg" alt="Py version">
  </a>
  <a href="https://hits.seeyoufarm.com">
    <img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftahiri-lab%2FiPhyloGeo_plus_plus&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false" alt="Hits">
  </a>
</p>

<h2 align="center">Multi-platform Application for Analyzing Phylogenetic Trees with Climatic Parameters</h2>

<table style="width:100%; border: none;">
  <tr>
    <td style="vertical-align: top; width: 50%;">
      <details open>
        <summary><strong>Table of Contents</strong></summary>
        <ol>
          <li><a href="#about-the-project">About the Project</a></li>
          <li><a href="#features">Features</a></li>
          <li><a href="#installation">Installation</a></li>
          <li><a href="#usage">Usage</a></li>
          <li><a href="#project-structure">Project Structure</a></li>
          <li><a href="#contributing">Contributing</a></li>
          <li><a href="#references">References</a></li>
          <li><a href="#contact">Contact</a></li>
        </ol>
      </details>
    </td>
    <td style="vertical-align: top; width: 50%; text-align: center;">
      <img src="/img/readme-pic/Main.png" alt="Main" style="max-width: 100%; border-radius: 8px;">
    </td>
  </tr>
</table>

# About the Project

The iPhyloGeo++ tool is an advanced bioinformatics application specifically designed for the integration and analysis of phylogeographic data. By leveraging both genetic and climatic information, it offers comprehensive insights into the evolutionary processes and geographical distribution of various species.

# Features

- **Cross-Platform Compatibility:** Compatible with Windows, macOS, and Linux.
- **Comprehensive Data Integration:** Merges genetic sequences with climatic data for robust analysis.
- **Intuitive User Interface:** Developed with PyQt6 to ensure ease of use.
- **Advanced Visualization Tools:** Provides visualization of phylogenetic trees and climatic data on interactive maps.
- **Enhanced Comparative Analysis:** Facilitates the comparison of different phylogenetic trees.

# [Installation](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/Getting-Started)

## 1. Clone the repository

```sh
git clone https://github.com/tahiri-lab/iPhyloGeo_plus_plus.git
cd iPhyloGeo_plus_plus
```

## 2. Install Poetry

Refer to the official installation guide [here](https://python-poetry.org/docs/#installation)

### Install with the [official installer](https://python-poetry.org/docs/#installing-with-the-official-installer)

#### Windows:

```ps
(Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | py -
```

#### Linux:

```sh
curl -sSL https://install.python-poetry.org | python3 -
```

### Install with `pipx`

```sh
pipx install poetry
```

## 2. Install the dependencies

Once `Poetry` installed, you can install the project's dependencies.
Go in the project directory and run the commands

```sh
# This will ensure that the virtual env will be created in the project directory
poetry config virtualenvs.in-project true

# Install the dependencies
poetry install
```

## 3. Run the Application

Once Poetry installed with the project's dependencies, you can activate the virtual env or run directly with poetry.

```sh
# Using the "vanilla" python way.
./.venv/Scripts/activate

# Using Poetry
iex (poetry env activate)

# After activating the environment
python ./scripts/main.py

# Alternatively, you can run the script directly without activating the venv.
poetry run python ./scripts/main.py
```

# [Usage](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/Discover-the-application)

## Loading Genetic Data

<p align="center"><img src="./img/readme-pic/Genetics.gif" alt="Genetic demonstration"></p>

1. **Navigate to File Browser on the Genetic Page:**

- Access the genetic data interface through the File Browser tab.
- Select and Load Your Fasta File:

To help you, you can go to the <a href="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/Discover-the-application">Wiki</a> to understand how to navigate into the iPhyloGeo application.

2. **Choose your Fasta file containing the genetic sequences. Supported formats should be specified (e.g., .fasta, .fa).**

- Ensure the file adheres to the correct format and structure.
- Perform Sequence Alignment, Statistics, and Generate Genetic Trees:

3. **Sequence Alignment:**

- Utilize built-in tools for aligning sequences, detailing available algorithms (e.g., MUSCLE, ClustalW).
- Statistics: Generate statistics such as nucleotide frequencies, sequence length distribution, and GC content.
- Genetic Trees: Construct phylogenetic trees using methods like Neighbor-Joining, Maximum Likelihood, or Bayesian inference. Visualize trees with options for customization (e.g., color-coding branches, annotating clades).

## Loading Climatic Data

<p align="center"><img src="./img/readme-pic/Climat.gif" alt="Climatic demonstration"></p>

1. **Navigate to File Browser on the Climatic Page:**

- Access the climatic data interface through the File Browser tab.
- Select and Load Your CSV File Containing Climatic Data:
- Choose your CSV file with climatic information. Supported data formats and required structure should be clarified.

2. **View the Generated Maps, Data Tables, Statistics, and Climatic Trees as Needed:**

- Display and interact with the visual representations of the climatic data, including maps, tables, and statistical summaries.

## Displaying Results

<p align="center"><img src="./img/readme-pic/Result.gif" alt="Results demonstration"></p>

1. **Navigate to the Results Page:**

- Access the results interface.

2. **Adjust the Parameters as Needed:**

- Modify settings to refine the analysis.

3. **Click on Submit to View the Phylogenetic Results:**

- Generate and display the results based on the input data and parameters.

4. **Navigate to the Stats Button for Phylogenetic Trees Visualization:**

- Use the stats button to visualize the phylogenetic trees and related statistics.

To help you, you can follow the <a href="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/Tutorial">Tutorial</a> part of the Wiki to achieve all these steps.

# Project Structure

This project is organized into several key directories to help you navigate and understand the codebase.

- **img/:** Contains images used by the README and the application.
- **datasets/:** Includes sample data for testing purposes.
- **scripts/:** Houses the Python files for the project.
- **poetry.lock:** List of dependencies.
- **scripts/main.py:** Main application entry point.

# Contributing

We welcome contributions to iPhyloGeo++. Please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Set up pre-commit hooks (`pre-commit install`).
4. Commit your changes (`git commit -m 'Add new feature'`).
5. Push to the branch (`git push origin feature-branch`).
6. Create a new Pull Request.

# ✔️ References

💡 If you are using our algorithm in your research, please cite our recent papers:

1️⃣ Li, W. & Tahiri, N. (2024). Host-Virus Cophylogeny Trajectories: Investigating Molecular Relationships between Coronaviruses and Bat Hosts. [Viruses, 16(7), p.1133](https://www.mdpi.com/1999-4915/16/7/1133)

2️⃣ Gagnon, J. & Tahiri, N. (2024). Ecological and Spatial Influences on the Genetics of Cumacea (Crustacea: Peracarida) in the Northern North Atlantic.
[Proceeding in SciPy 2024, Tacoma, WA, USA](https://proceedings.scipy.org/articles/NVYF1037)

3️⃣ Li, W. & Tahiri, N. (2023). aPhyloGeo-Covid: A Web Interface for Reproducible Phylogeographic Analysis of SARS-CoV-2 Variation using Neo4j and Snakemake.
[Proceeding in SciPy 2023, Auxtin, TX, USA](https://conference.scipy.org/proceedings/scipy2023/pdfs/nadia_tahiri.pdf)

:four: Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.
[Proceeding in SciPy 2022, Auxtin, TX, USA](https://conference.scipy.org/proceedings/scipy2022/pdfs/nadia_tahiri.pdf)

# Contact

Please email us at: <Nadia.Tahiri@USherbrooke.ca> for any questions or feedback.

```

```

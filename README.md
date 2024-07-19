﻿<h1  align="center"> 🌳 iPhyloGeo++ <p align='center'> 
        [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 
        [![Contributions](https://img.shields.io/badge/contributions-welcome-blue.svg)](https://pysd.readthedocs.io/en/latest/development/development_index.html)
        [![Py version](https://img.shields.io/pypi/pyversions/pysd.svg)](https://pypi.python.org/pypi/pysd/)
        [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftahiri-lab%2FaPhyloGeo_plus_plus&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
        </p>


<h2  align="center">Multi-platform application for analyze phylogenetic trees with climatic parameters</h2>

<table>
<tr>
<th align="left">
<img width="441" height="1">
<p> 
<small>
<details open>
  <summary>Table of Contents</summary>
  <ol>
    <li>
    About the project
    </li>
    <li>Features
    </li>
    <li>Installation</li>
    <li>Usage</li>
        <li>Project Structure</li>
   <li>Contributing</li>
   <li>References</li>
   <li>Contact</li>

  </ol>
</details>
</small>
</p>
</th>
<th align="center">
<img width="441" height="1">
<p> 
<small>
<img src="/img/readme-pic/Main.png"  alt="1">
</small>
</p>
</th>
</tr>
</table>


# About the project

`iPhyloGeo++` is a sophisticated bioinformatics tool designed for integrating and analyzing phylogeographic data. This application leverages genetic and climatic data to provide comprehensive insights into the evolutionary and geographical distribution of species.

# Features
- **Cross-Platform Support:** Available on Windows, macOS, and Linux.
- **Data Integration:** Combines genetic sequences with climatic data.
- **User-Friendly Interface:** Built with PyQt5 for an intuitive experience.
- **Visualization Tools:** Displays phylogenetic trees and climatic data on maps.
- **Comparative Analysis:** Enables comparison of phylogenetic trees.

# Installation
**1- Clone the repository**
```sh
git clone https://github.com/tahiri-lab/iPhyloGeo_plus_plus.git
cd iPhyloGeo_plus_plus
```

**2- Set Up a Virtual Environment**
```sh
python3 -m venv iPhyloGeo_env
source iPhyloGeo_env/bin/activate  # On Windows use `iPhyloGeo_env\Scripts\activate`
```

**3- Install Dependencies:**
```sh
pip install -r requirements.txt
```

**4-Run the Application:**
```sh
python3 main.py
```

# Usage
## Loading Genetic Data
<p align="center">
  <img src="./img/other/genetic.gif" alt="Genetic demonstration">
</p>

1. Navigate to *File Browser* in the genetic page.
2. Select and load your Fasta file.
3. Perform sequence alignment, statistics and genetic trees as needed.
   
## Loading Climatic Data
<p align="center">
  <img src="./img/other/climatic.gif" alt="Climatic demonstration">
</p>

1. Navigate to *File Browser* in the climatic page.
2. Select and load your CSV file containing climatic data.
3. View the generated maps, corresponding data tables, statistics and climatic trees as needed.
   
## Displaying Results
<p align="center">
  <img src="./img/other/results.gif" alt="Results demonstration">
</p>

1. Navigate to the results page.
2. Adjust the parameters as needed.
3. Click on *Submit* to view the phylogenetic results
4. navigate to the *stats* button for the phylogenetic trees visualization.

# Project Structure
This project is organized into several key directories to help you navigate and understand the codebase.

- **img/:** Contains images used by the README and the application.
- **datasets/:** Includes sample data for testing purposes.
- **scripts/:** Houses the Python files for the project.
- **requirements.txt:** List of dependencies
- **scripts/main.py:** Main application entry point


# Contributing
We welcome contributions to iPhyloGeo++. Please follow these steps:
1. Fork the repository.
2. Create a new branch ( ```sh git checkout -b feature-branch``` ).
3. Commit your changes ( ```sh git commit -am 'Add new feature' ``` ).
4. Push to the branch ( ```sh git push origin feature-branch``` ).
5. Create a new Pull Request.



# References

1️⃣ Calculation of distance between phylogenetic tree: **Least Square metric**
+ [Cavalli-Sforza, L. L., & Edwards, A. W. (1967). Phylogenetic analysis. Models and estimation procedures. American journal of human genetics, 19(3 Pt 1), 233.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/)
+ [Felsenstein, J. (1997). An alternating least squares approach to inferring phylogenies from pairwise distances. Systematic biology, 46(1), 101-111.](https://pubmed.ncbi.nlm.nih.gov/11975348/)
+ [Makarenkov, V., & Lapointe, F. J. (2004). A weighted least-squares approach for inferring phylogenies from incomplete distance matrices. Bioinformatics, 20(13), 2113-2121.](https://pubmed.ncbi.nlm.nih.gov/15059836/)

2️⃣ Calculation of distance between phylogenetic tree: **Robinson-Foulds metric**
+ [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)
    
3️⃣ Dataset full description: **Analysis of genetic and climatic data of SARS-CoV-2**
+ [Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.](https://conference.scipy.org/proceedings/scipy2022/nadia_tahiri.html)

# 📧 Contact
Please email us at : <Nadia.Tahiri@USherbrooke.ca> for any question or feedback.


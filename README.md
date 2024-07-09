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
      <a href="#about-the-project">About the project</a>
    </li>
    <li>
      <a href="#Installation">Installation</a>
      <ul>
        <li><a href="#Linux-UNIX-and-Mac-OS-versions">Linux/UNIX and Mac OS versions</a></li>
      </ul>
    </li>
    <li> How to use</li>
      <ul>
        <li><a href="#Getting-genetic-data">Getting genetic data</a></li>
        <li><a href="#Getting-climatic-data">Getting climatic data</a></li>
        <li><a href="#Display-results">Display results</a></li>
      </ul>
    <li>
      <a href="#contact">Contact</a>
    </li>
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


# 📝 About the project

`iPhyloGeo++` is a bioinformatics pipeline dedicated to the analysis of phylogeography. `iPhyloGeo++` is an open-source multi-platform application designed by the team of Professor [Nadia Tahiri](https://tahirinadia.github.io/) (University of Sherbrooke, Quebec, Canada). It is implemented in Python. This tool can be used to obtain trees from climatic data of the regions where the samples have been collected. Those climatic trees are then used for topological and evolutionary comparison against phylogenetic trees from multiple sequence alignments (MSAs) using the [Least Square (LS) metric](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/). MSAs that yield trees with a significant `LS` value are then optionnally saved in folders with their respective tree. The `output.csv` file contains the informations of all the significant MSAs informations (see Worflow Section for more details).


# ⚒️ Installation

## Linux UNIX and Mac OS versions
`iPhyloGeo++` is available as a Python script.

### Prerequisites
💡 Before using this program, make sure that you have installed all the necessary libraries for it to work properly. To do this, simply type the following command:

```
pip3 install -r requirements.txt
```

### Python script
A `requirements.txt` file containing all required libraries is available in the GitHub repository.

⚠️ Assuming Python 3.8 or higher is installed on the machine, the script should run well with the libraries installed.

<u>Here is an example of how to run the script in Linux/UNIX or Mac OS:</u>
1. After downloading the source code, go to the folder containing `main.py`.
2. If you do not have `virtualenv` installed, run `python3 -m pip install --user virtualenv`
3. Create a new virtual environment (venv) in your terminal using `python3 -m venv iPhyloGeo++_env`.
4. Still in the terminal, enter the new venv using `source iPhyloGeo++_env/bin/activate`.
5. Install the required libraries using `pip install -r requirements.txt`.
6. Launch iPhyloGeo using `python3 main.py`.

<!--You can also launch the package using the `make` command from your terminal when you are in the `root`. This command will use the `Makefile` to run the script. If you use the command `make clean`, it will erase the `output.csv` file previously created with the first command.-->


# Getting started

## Menu 
When you first launch the application you can see a window with multiple buttons at the top and some informations about how to use `iPhyloGeo++`.

![group_page_step](/img/readme-pic/Iphy-Menu.png)

The menu is composed of 6 buttons : 

### Home button : 
![group_page_step](/img/readme-pic/Home-Button.png)
### Genetic data button : 
![group_page_step](/img/readme-pic/Genetic-Button.png)
### Climatic data button : 
![group_page_step](/img/readme-pic/Climatic-Button.png)
### Result button : 
![group_page_step](/img/readme-pic/Result-Button.png)
### Help button : 
![group_page_step](/img/readme-pic/Help-Button.png)
### Light/Dark button : 
![group_page_step](/img/readme-pic/LightDark-Button.png)

## Genetic Data
In the genetic data tab, the user can load a sequence file (only Fasta files are accepted) and alignment will be done when using the `Start` button in the `sequence alignment` tab.


### A typical workflow would look like this:

![group_page_step](/img/readme-pic/Aphy-Use.png)

1. Select your sequence file in an accepted format (Fasta) through the **File browser** button.
2. Once the file is loaded, the sequence will appear in the main window.
3. Go to the sequence alignment page through the **Sequence Alignment** button. 
4. select the alignment method that you want to use and press the **Start button**.
  - The sequence alignment process might be CPU demanding so a good PC configuration can be necessary, with an Mac M1 chip, the whole process can take 2 to 3 minutes.
  - When the process is done the sequence aligned will be displayed on the window

Below is a summary of the steps presented:

![group page gif](/img/readme-pic/geneticData-gif.gif)

## Climatic Data
After selecting the required file for the climatic data section, a tab and a map will automatically be generated.

1. Select your climatic file in an accepted format (csv) through the **File browser** button.
2. Once the file is loaded, the tab and the map will appear in the main window.

Below is a summary of the steps presented:

![group page gif](/img/readme-pic/climaticData.gif)

## Display results
With both genetic and climatic data we can now move forward to the last section, the **Results** tab.

This feature will display the results through the **submit** button:

![group page gif](/img/readme-pic/results1.gif)

# ✔️ References

1️⃣ Calculation of distance between phylogenetic tree: `Least Square metric`
+ [Cavalli-Sforza, L. L., & Edwards, A. W. (1967). Phylogenetic analysis. Models and estimation procedures. American journal of human genetics, 19(3 Pt 1), 233.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/)
+ [Felsenstein, J. (1997). An alternating least squares approach to inferring phylogenies from pairwise distances. Systematic biology, 46(1), 101-111.](https://pubmed.ncbi.nlm.nih.gov/11975348/)
+ [Makarenkov, V., & Lapointe, F. J. (2004). A weighted least-squares approach for inferring phylogenies from incomplete distance matrices. Bioinformatics, 20(13), 2113-2121.](https://pubmed.ncbi.nlm.nih.gov/15059836/)

2️⃣ Calculation of distance between phylogenetic tree: `Robinson-Foulds metric`
+ [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)
    
3️⃣ Dataset full description: `Analysis of genetic and climatic data of SARS-CoV-2`
+ [Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.](https://conference.scipy.org/proceedings/scipy2022/nadia_tahiri.html)

# 📧 Contact
Please email us at : <Nadia.Tahiri@USherbrooke.ca> for any question or feedback.


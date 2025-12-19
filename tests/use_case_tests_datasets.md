# Use Case Testing: Datasets

This set of tests can be used to confirm that iPhyloGeo++ works correctly with a greater variety of datasets.

## Setup

### Install or Update iPhyloGeo++

Make sure you have the latest version of iPhyloGeo++ installed. Follow the appropriate install guide if necessary: [Windows](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/3.2-Installation-(Windows)) [macOS](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/3.1-Installation-(MAC)) [Linux](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/3.3-Installation-(Linux))

### Get the Datasets

Download the following files:
1. [Cumacea.fasta](https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/Cumacea/Cumacea.fasta)
2. [Cumacea.csv](https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/Cumacea/Cumacea.csv)
3. [The_37_climate.csv](https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/covid/The_37_climate.csv)
4. [The_37_seq.fasta](https://github.com/tahiri-lab/aPhyloGeo/blob/main/datasets/covid/The_37_seq.fasta)

If you installed the entire project by cloning the repository or copying the entire source code, you already have **geo_with_loc.csv**, **simplot.csv**, **simplot.fasta** and **small_seq.fasta** They are located in the **datasets** folder.

If you only installed the Windows executable, download them from [this folder](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/datasets/) 

## Tests

### **1.1 Genetic Analysis**

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Genetic-Button.png)

**Action:** Click the **Genetic** button.  
**Expected Result:** The Genetic **Get Started** screen is displayed. On the left, the **File Browser** and **Clear** buttons are active while others are greyed out.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/genetic_screen.png" alt="Home screen" width="500"/>

**Action:** Click the **File Browser** button.
**Expected Result:** A window opens allowing you to select a FASTA file to upload the genetic data.

**Action:** Browse to and select **Cumacea.fasta**.
**Expected Result:** The file selection window closes and the contents of the file are displayed. Generic text is displayed in green and genetic sequences are displayed in green (A), blue (C), red (G) and black (T). A small scrollbar appears on the right.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/fasta.png" alt="Fasta File displayed within iPhyloGeo++" width="500"/>

**Action:** Click the **Sequence Alignment** tab.
**Expected Result:** A blank screen with the title "Genetic Sequence" and two green buttons, a "Play" button and a "Settings" (gear icon) button, is displayed.

**Action:** Click the green **Settings** (gear icon) button.
**Expected Result:** A popup window titled "Parameter Dialog" appears.

**Action:** Select the following settings:
* Bootstrap Threshold: **0**
* Step Size: **100**
* Bootstrap Amount: **100**
* Alignment Method: **CLUSTALW**
* Fit Method: **Narrow-fit prevent elongation with gap when possible**
* Tree Type: **FastTree application**
* Rate Similarity: **90**
* Method Similarity: **Jaro similarity**
**Expected Result:** The user is able to select the above settings.

**Action:** Click the **Save Settings** button.
**Expected Result:** The popup window disappears.

**Action**: Click the **Play** button.
**Expected Result:** A console window and a loading dialog box may appear. Then, the genetic sequences are displayed as in the image below. Note that the first colomn of the sequence (to the immediate right of the species labels) should contain only instances of T (black letter T on a red background) and the last column of the sequence should contain both T (red background) and C (yellow background)

TODO complete the test and add screenshots of the results once [issue 30](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/issues/30) is fixed.

## 2. External Documentation

Additional test references are available on GitHub:

- **Application Discovery:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/5.-Discover-the-application  

- **Tutorial:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/6.-Tutorial
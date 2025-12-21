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

## 1 Tests

### **1.1 Genetic Analysis**

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Genetic-Button.png)

**Action:** Click the **Genetic** button.  
**Expected Result:** The Genetic **Get Started** screen is displayed. On the left, the **File Browser** and **Clear** buttons are active while others are greyed out.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/genetic_screen.png" alt="Genetic Get Started screen" width="500"/>

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
* Alignment Method: **MUSCLE**
* Fit Method: **Narrow-fit prevent elongation with gap when possible**
* Tree Type: **FastTree application**
* Rate Similarity: **90**
* Method Similarity: **Jaro similarity**
**Expected Result:** The user is able to select the above settings.

**Action:** Click the **Save Settings** button.
**Expected Result:** The popup window disappears.

**Action:** Click the **Play** button.
**Expected Result:** A console window and a loading dialog box may appear. Then, the genetic sequences are displayed.

**Action:** Set the **Window size** to 35 and the **Starting position** to 1.
**Expected Results:** If the settings have been changed, the display will update. The result should match the image below:

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/cumacea_alignment.png" alt="Genetic sequence alignment for Cumacea" width="500"/>

### **1.2 Climatic Analysis**

#### Test 1.2.1

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Climatic-Button.png)

**Action:** Click the **Climate** button.  
**Expected Result:** The Climatic **Get Started** screen is displayed. On the left, the **File Browser** and **Clear** buttons are active while others are greyed out.

**Action:** Set the Max Correlation to 1,0000 and the Min Variance to 0,0000.
**Expected Result:** You are able to set the Max Correlation and Min Variance values.

**Action:** Click the **File Browser** button.  
**Expected Result:** A window opens allowing you to select a CSV file to upload the climatic data.

**Action:** Browse to and select **cumacea.csv**.
**Expected Result:** The contents of the file are displayed in the form of a two-dimensional table with 8 columns, with a map underneath, as in the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/cumacea_table_map.png" alt="A table and a map representing climatic data" width="500"/>

**Action:** Click the **Statistics** button on the left panel.
**Expected Result:** On the left, the title "Statistics" is displayed, with the subtitle "Generate your graph", three dropdown menus, "Insert X axis data", "Insert Y axis data" and "Choose plot type", and a green save button. On the right, a graph is displayed (it may take several seconds to load)

**Action:** Select Correlation as the plot type.
**Expected Result:** The Insert X axis data and Insert Y axis data dropdown menus are disabled. The graph displayed matches the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/correlation.png" alt="A correlation map" width="500"/>

**Action:** Click the **Climatic Tree** button on the top menu or the left panel.  
**Expected Result:** A climatic tree is displayed (it may take a few seconds), with a Gear/Settings button at the top left and a Save button at the top right. A "Variable to plot" dropdown menu is displayed to the right of the gear icon.

**Action:** Click the **Gear/Settings** button.
**Expected Result:** A popup window titled Preferences and matching the screenshot below is displayed.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/climatic_tree_prefs.png" alt="A preferences window" width="150"/>

**Action:** Select "black" as the Label color, "blue" as the Edge color, "Hierarchical Horizontal" as the Layout Option, and "Tree View" as the View type. Of the 4 checkboxes, check "Use leaf names" and leave the three others unchecked. Click **Save**.
**Expected Result:** The Preferences window disappears. The tree may be updated.

**Action:** Select "O2-saturation_ground" as the **Variable to plot**.
**Expected Result:** The climatic tree is updated to match the screenshot below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/o2tree.png" alt="A climatic tree" width="500"/>

#### Test 1.2.2

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Climatic-Button.png)

**Action:** Click the **Climate** button.  
**Expected Result:** The Climatic **Get Started** screen is displayed. On the left, the **File Browser** and **Clear** buttons are active while others are greyed out.

**Action:** Set the Max Correlation to 0,8000 and the Min Variance to 0,0001.
**Expected Result:** You are able to set the Max Correlation and Min Variance values.

**Action:** Click the **File Browser** button.  
**Expected Result:** A window opens allowing you to select a CSV file to upload the climatic data.

**Action:** Browse to and select **The_37_climate.csv**.
**Expected Result:** The contents of the file are displayed in the form of a two-dimensional table with 6 columns, with a map underneath, as in the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/37climatic_data.png" alt="A table and a map representing climatic data" width="500"/>

**Action:** Click the **Statistics** button on the left panel.
**Expected Result:** On the left, the title "Statistics" is displayed, with the subtitle "Generate your graph", three dropdown menus, "Insert X axis data", "Insert Y axis data" and "Choose plot type", and a green save button. On the right, a graph is displayed (it may take several seconds to load)

**Action:** Select PRECTOTCORR for the X axis, ALLSKY_SFC_SW_DWN for the X axis and Scatter Plot as the plot type.
**Expected Result:** The graph displayed matches the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/37climatic_stats.png" alt="A scatter plot from the Climatic section" width="500"/>

**Action:** Click the **Climatic Tree** button on the top menu or the left panel.  
**Expected Result:** A climatic tree is displayed (it may take a few seconds), with a Gear/Settings button at the top left and a Save button at the top right. A "Variable to plot" dropdown menu is displayed to the right of the gear icon.

**Action:** Click the **Gear/Settings** button.
**Expected Result:** A popup window titled Preferences is displayed.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/climatic_tree_prefs.png" alt="A preferences window" width="150"/>

**Action:** Select "black" as the Label color, "red" as the Edge color, "Radial" as the Layout Option, and "Tree View" as the View type. Check the three first checkboxes, leaving Show branch lengths unchecked. Click **Save**.
**Expected Result:** The Preferences window disappears. The tree may be updated.

**Action:** Select "T2M" as the **Variable to plot**.
**Expected Result:** The climatic tree is updated to match the screenshot below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/37climatic_tree.png" alt="A climatic tree" width="500"/>

### **1.3 Results Section**

**Setup:**
- Click the **Genetic** button, then the **File Browser** button. Provide **small_seq.fasta**.
- Click the **Alignment**, button, then the **Play** button. Wait for the sequences to be displayed.
- Click the **Climatic** button, then the **File Browser** button. Provide **geo_with_loc.csv**.

**Test:**

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Result-Button.png)

**Action:** Click the **Results** button (immediately to the left of the Help button).  
**Expected Result:** A blank screen with the title "Results" is displayed.

**Action:** Click the **Settings** button.
**Expected Result:** A dialog box titled Parameters appears, like on the image below. It contains three dropdown menus (Calculus method, Bootstrap threshold and Metric threshold). At the bottom are three buttons: Reset, OK and Cancel.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/results_params.png" alt="Parameters dialog box" width="500"/>

**Action:** Set the **Calculus method** to "Robinson and Foulds", the **Bootstrap threshold** to 0 and the **Metric threshold** to 50, then click **OK**.
**Expected Result:** The dialog box disappears.

**Action:** Click the **Statistics** button on the left panel.
**Expected Result:** On the left are displayed the title "Statistics" and a phylogenetic tree. On the right are displayed a green Download button, dropdown menus titled "Condition" and "trees" and a graph for the selected condition.

**Action:** Click the **Map** tab.
**Expected Result:** A loading bar may appear. Then, a graph with a phylogenetic tree on the left and species coordinates on the right appears, matching the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/geoeuclidmap.png" alt="Map screen" width="500"/>

**Action:** Click the **Save** button.
**Expected Result:** A screen with a table is displayed. The results are saved as scripts/results/output.csv

**Action:** Open scripts/results/output.csv using spreadsheet software.
**Expected Result:** The file’s contents should match the table in iPhyloGeo++.

### **1.4 Full Pipeline**

Tests 1.4.1, 1.4.2 and 1.4.3 should be performed in succession as 1.4.1 and 1.4.2 are the setup for 1.4.3.

#### 1.4.1 Climatic Section

TODO

#### 1.4.2 Genetic Section

TODO

TODO simplot

## 2 External Documentation

Additional test references are available on GitHub:

- **Application Discovery:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/5.-Discover-the-application  

- **Tutorial:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/6.-Tutorial
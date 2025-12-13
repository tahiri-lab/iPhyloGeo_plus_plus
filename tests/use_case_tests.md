# Use Case Testing

## Setup

### Install or Update iPhyloGeo++

Make sure you have the latest version of iPhyloGeo++ installed. Follow the appropriate install guide if necessary: [Windows](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/3.2-Installation-(Windows)) [macOS](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/3.1-Installation-(MAC)) [Linux](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/3.3-Installation-(Linux)

### Get the Dataset

If you installed the entire project by cloning the repository or copying the entire source code, you already have **simplot.csv** and **simplot.fasta**. They are located in the **datasets** folder.

If you only installed the Windows executable, download them from GitHub: [csv](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/datasets/simplot.csv) [fasta](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/datasets/simplot.fasta)

## 1. Use Cases

### **1.1. Home Button**
**Action:** Click the **Climatic** button (thermometer icon).
**Expected Result:** The home screen is replaced by the Climatic **Get Started** screen.

**Action:** Click the **Home** button.  
**Expected Result:** The home screen loads successfully. The title "Welcome to iPhyloGeo", the subtitle "Thank you for downloading our software", the main text "Here is your guide..." and the copyright mentions are displayed.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/home_page.png" alt="Home screen" width="500"/>

---

### **1.2. Genetic Analysis**

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Genetic-Button.png)

**Action:** Click the **Genetic** button.  
**Expected Result:** The Genetic **Get Started** screen is displayed. On the left, the **File Browser** and **Clear** buttons are active while others are greyed out.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/genetic_screen.png" alt="Home screen" width="500"/>

**Action:** Click the **File Browser** button.
**Expected Result:** A window opens allowing you to select a FASTA file to upload the genetic data.

**Action:** Browse to and select **simplot.fasta**.
**Expected Result:** The file selection window closes and the contents of the file are displayed. Generic text is displayed in green and genetic sequences are displayed in green (A), blue (C), red (G) and black (T). A small scrollbar appears on the right.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/fasta.png" alt="Fasta File displayed within iPhyloGeo++" width="500"/>

**Action:** Click the scrollbar on the right and bring it to the bottom.
**Expected Result:** The text scrolls all the way down. The last/bottom genetic sequence displayed is identified as "HPU63231 Helicobacter pylori".

**Action:** Click the **Sequence Alignment** tab.
**Expected Result:** A blank screen with the title "Genetic Sequence" and two green buttons, a "Play" button and a "Settings" (gear icon) button, is displayed.

**Action:** Click the **Fasta file** tab.
**Expected Result:** The genetic sequences are displayed once again. The scrollbar on the right is still at the bottom.

**Action:** Click the **Alignment** button on the left.
**Expected Result:** A blank screen with the title "Genetic Sequence" and two green buttons, a "Play" button and a "Settings" (gear icon) button, is displayed.

**Action:** Click the green **Settings** (gear icon) button.
**Expected Result:** A popup window titled "Parameter Dialog" appears.

**Action:** Select the following settings:
* Bootstrap Threshold: **0**
* Window Size: **35**
* Step Size: 100
* Bootstrap Amount: 100
* Alignment Method: MUSCLE
* Fit Method: WiderFit by elongating with Gap (starAlignment)
* Tree Type: FastTree application
* Rate Similarity: 90
* Method Similarity: Hamming distance
**Expected Result:** The user is able to select the above settings.

**Action:** Click the **Save Settings** button.
**Expected Result:** The popup window disappears.

**Action**: Click the **Play** button.
**Expected Result:** A console window and a loading dialog box may appear. Then, the genetic sequences are displayed as in the image below. Note that the first colomn of the sequence (to the immediate right of the species labels) should contain only instances of T (black letter T on a red background) and the last column of the sequence should contain both T (red background) and C (yellow background)

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/sequencepos1.png" alt="A sequence alignment for Helicobacter pylori, starting from position 1" width="500"/>

**Action**: At the bottom, enter 50 as the new window size: select the **Window size** selector, delete "35", input "50" and press the return key.
**Expected Result:** More columns are displayed, with the first and last column only containing T (red background).

**Action**: At the bottom, press the up arrow on the **Starting position** selector.
**Expected Result:** The columns have changed positions. The first and last displayed column now contain C (yellow background).

**Action**: Click the **Statistics** button on the left.
**Expected Result**: The title Alignment Chart is displayed, with a dropdown menu titled "Reference", a green Download button and a graph titled "Sequence Similarity Plot" matching the image below. It may take a few seconds for the dropdown menu to populate and the graph to appear.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/seqsimplot.png" alt="A sequence similarity plot for Helicobacter pylori" width="500"/>

**Action:** Hover over coordinates (10, 1) on the graph.
**Expected Result:** Hover displays of multiple different colors will appear (their number will depend on display size). The one at the top will read "Sequence: HPU63239 Helicobacter pylori" on the first line and "Similarity: 1.00" on the second line. On the horizontal axis, the value 10 will be highlighted.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/seqsimplot10-1.png" alt="A sequence similarity plot for Helicobacter pylori" width="500"/>

**Action:** Select another value in the **Reference** dropdown menu.
**Expected Result:** The graph will be updated.

**Action:** Click the **Download** button.
**Expected Result:** A window appears to allow the user to select a directory and file name to save the graph.

**Action:** Provide a directory and PNG file name of your choice. Save the file.
**Expected Result:** The graph is saved in the location you provided. It matches the colors and values of the graph displayed in iPhyloGeo++ (note that it may not have the same dimensions).

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/seqsimplot_comparison.png" alt="A sequence similarity plot for Helicobacter pylori" width="700"/>

**Action:** Click the **Genetic Tree** button.  
**Expected Result:** A console window may appear and disappear. Then, the title "Phylogenetic Tree", a Reference dropdown menu, a green Download button and a phylogenetic tree matching the image below should be displayed.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/0nt19nt.png" alt="A phylogenetic tree for Helicobacter pylori" width="500"/>

**Action:** In the **Reference** dropdown menu, select "40 nt 59 nt". 
**Expected Result:** The displayed phylogenetic tree will be updated to match the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/40nt59nt.png" alt="A phylogenetic tree for Helicobacter pylori" width="500"/>

---

### **1.3. Climate Button**  
![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Climatic-Button.png)

**Action:** Click the **Climate** button.  
**Expected Result:** The Climate module opens and all available options are displayed.

**Action:** Click the **Get Started** button.  
**Expected Result:** You are able to change values in the dropdown list.

**Action:** Click the **File Browser** button.  
**Expected Result:** A window opens allowing you to select a CSV file to upload the climatic data.

**Action:** Click the **Statistics** button on the left panel.  
**Expected Result:** Climate-related statistics are displayed.

**Action:** Click the **Climatic Data** button on the top menu.  
**Expected Result:** All climatic data are displayed.

**Action:** Select an option from the **dropdown list**.  
**Expected Result:** The corresponding data are displayed successfully.

**Action:** Click the **Climatic Tree** button on the top menu or the left panel.  
**Expected Result:** The climatic tree is generated and displayed.

**Action:** Click the **Download** button.  
**Expected Result:** A download window opens successfully.

**Action:** Click the **dropdown list** under the Climatic Tree section.  
**Expected Result:** The corresponding data load correctly.

**Action:** Click the **Settings** button under the Climatic Tree section.  
**Expected Result:** All available settings options are displayed for modification.

---

### **1.4. Results Section**

#### Setup

Click the **Genetic** button, then the **File Browser** button. Provide **simplot.fasta**.

Click the **Alignment**, button, then the **Play** button. Wait for the sequences to be displayed.

Click the **Climatic** button, then the **File Browser** button. Provide **simplot.csv**.

#### Test

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Result-Button.png)

**Action:** Click the **Results** button (immediately to the left of the Help button).  
**Expected Result:** A blank screen with the title "Results" is displayed.

**Action:** Click the **Settings** button.
**Expected Result:** A dialog box titled Parameters appears, like on the image below. It contains three dropdown menus (Calculus method, Bootstrap threshold and Metric threshold). At the bottom are three buttons: Reset, OK and Cancel.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/results_params.png" alt="Parameters dialog box" width="500"/>

**Action:** Set the **Calculus method** to "Robinson and Foulds", the **Bootstrap threshold** to 5 and the **Metric threshold** to 90, then click **OK**.
**Expected Result:** The dialog box disappears.

**Action:** Click the **Settings** button.
**Expected Result:** The dialog box reappears, with values "Robinson and Foulds", 5 and 90 displayed.

**Action:** Click the **Reset** button.
**Expected Result:** The selected values remain "Robinson and Foulds", 5 and 90.

**Action**: Set the **Bootstrap threshold** to 10, then click the **Reset** button.
**Expected Result:** The Boostrap threshold is set back to 5.

**Action**: Set the **Bootstrap threshold** to 10, then click the **Cancel** button.
**Expected Result:** The dialog box disappears.

**Action:** Click the **Settings** button.
**Expected Result:** The dialog box reappears, with values "Robinson and Foulds", 5 and 90 displayed.

**Action:** Set the **Calculus method** to "Least square", the **Bootstrap threshold** to 180 and the **Metric threshold** to 60, then click **OK**.
**Expected Result:** The dialog box disappears.

**Action:** Click the **Statistics** tab (between Results and Map).
**Expected Result:** On the left are displayed the title "Statistics" and a phylogenetic tree. On the right are displayed a green Download button, dropdown menus titled "Condition" and "trees" and a graph for the selected condition (in this case, the graph will be empty), matching the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/statistics.png" alt="Statistics screen" width="500"/>

**Action:** Set "ALLSKY_SFC_SW_DWN" as the Condition and "40 nt 59 nt" as the trees value.
**Expected Result:** The graph on the left will update, matching the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/statistics40nt59nt.png" alt="Statistics screen" width="500"/>

**Action:** Click the **Map** tab.
**Expected Result:** A loading bar may appear. Then, a graph with a phylogenetic tree on the left and species coordinates on the right appears, matching the image below.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/map.png" alt="Map screen" width="500"/>

**Action:** Click the **Save** button.
**Expected Result:** A screen with a list of headers is displayed, matching the image below. The results are saved in the project’s scrips/results directory.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/save.png" alt="Save screen" width="500"/>

TODO specify how to verify the contents of the results folder.

---

### **1.5. Help Button**  
![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Help-Button.png)

**Action:** Click the **Help** button.  
**Expected Result:** The help/documentation window opens correctly.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/help_popup.png" alt="Help window" width="500"/>

---

### **1.6. Light Mode Button**
**Action:** Toggle the **Light Mode** button. (moon/sun button)
**Expected Result:** The application switches successfully between light and dark themes.
 
<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/darkmode.png" alt="Dark mode" width="200"/> <img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/lightmode.png" alt="Light mode" width="200"/>
---


---

## 2. External Documentation

Additional test references are available on GitHub:

- **Application Discovery:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/5.-Discover-the-application  

- **Tutorial:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/6.-Tutorial

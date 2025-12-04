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

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/home_page.png" alt="Home screen" width="350"/>

---

### **1.2. Genetic Button**

![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Genetic-Button.png)

**Action:** Click the **Genetic** button.  
**Expected Result:** The Genetic **Get Started** screen is displayed. On the left, the **File Browser** and **Clear** buttons are active while others are greyed out.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/genetic_screen.png" alt="Home screen" width="350"/>

**Action:** Click the **File Browser** button.
**Expected Result:** A window opens allowing you to select a FASTA file to upload the genetic data.

**Action:** brows to this location **iPhyloGeo_plus_plus/datasets** and select  **Fasta File**.

**Action:** Click the **setting** button to run the life. 
**Action:** you can use the default seetings or change them **(ex.: "Select the following settings: * Bootstrap Threshold: 0; Window Size: 20 (...)**.
**Expected Result:** see video down below.
![Genetics](https://github.com/user-attachments/assets/6187329c-506b-4a9a-8285-7c83143ab751)

**Action:** Click the **Sequence Alignment** button on the top menu or **Alignment** on the left panel.  
**Expected Result:** The *Play Sequence* button and the *Settings* button appear.

**Action:** Click the **Play Sequence** button.  
**Expected Result:** The *Statistics* and *Genetic Tree* options become unlocked and available.

**Action:** Click the **Settings** button.  
**Expected Result:** All available settings options are displayed for modification.

**Action:** Click the **Statistics** button.  
**Expected Result:**  
- Statistical outputs are generated and displayed.  
- Selecting items from the dropdown list displays the correct data.  
- Clicking the **Download** button opens the download window successfully.

**Action:** Click the **Genetic Tree** button.  
**Expected Result:**  
- The genetic tree is generated and displayed successfully.  
- Selecting items from the dropdown list updates the displayed tree or data.  
- Clicking the **Download** button opens the download window successfully.

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

### **1.4. Results Button**  
![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Result-Button.png)

**Action:** Click the **Results** button.  
**Expected Result:**  
- The Results view is displayed.  
- This menu **only** displays results (no uploads or processing actions).

**Action:** Click the **Results** button once both datasets have been uploaded.  
**Expected Result:**  
- The Results section becomes accessible only when both genetic and climatic data are available.

---

### **1.5. Help Button**  
![image](https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/readme-pic/Help-Button.png)

**Action:** Click the **Help** button.  
**Expected Result:** The help/documentation window opens correctly.

<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/help_popup.png" alt="Help window" width="350"/>

---

### **1.6. Light Mode Button**
**Action:** Toggle the **Light Mode** button. (moon/sun button)
**Expected Result:** The application switches successfully between light and dark themes.
 
<img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/darkmode.png" alt="Dark mode" width="200"/> <img src="https://github.com/tahiri-lab/iPhyloGeo_plus_plus/blob/main/img/usecase-pic/lightmode.png" alt="Light mode" width="200"/>
---

### **1.7. Clear Button**
**Action:** Click the **Clear** button.  
**Expected Result:** All displayed data and information are deleted/reset.


---

## 2. External Documentation

Additional test references are available on GitHub:

- **Application Discovery:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/5.-Discover-the-application  

- **Tutorial:**  
  https://github.com/tahiri-lab/iPhyloGeo_plus_plus/wiki/6.-Tutorial

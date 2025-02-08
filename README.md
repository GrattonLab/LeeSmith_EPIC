# 1. EPIC (Extended Precision measurement of Inhibitory Control)

### 1.1 Description
- This project examines how precise individual-level estimates of inhibitory control measures can be when utilizing extended within-subject data (exceeding 5,000 trials per participant per task). The dataset is publicly available at [OSF](https://osf.io/jk9nb/?view_only=bfc4c56718334581b267ad2e6f970d74).

### 1.2 Citation
- The codes provided here are based on the methods described in the following article: <br/>
Lee, H. J., Smith, D. M., Hauenstein, C., Dworetsky, A., Kraus, B. T., Dorn, M., Nee, D. E., & Gratton, C. (2025). [Precise individual measures of inhibitory control](https://osf.io/preprints/psyarxiv/rj2bu_v1). *Nature Human Behaviour*

### 1.3 Features
- All analyses can be replicated with the provided datasets and scripts.

---

# 2. System Requirements
All MATLAB scripts are written and tested in MATLAB 2021b. 

- **MATLAB System Requirements:** [View Requirements](https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
- **MATLAB Installation:** [Install MATLAB](https://www.mathworks.com/help/install/install-products.html)

*Note:* The installation can take up to 30 minutes, depending on your computer model.

**Required for Experimental Tasks:**
- **Psychtoolbox** for MATLAB: [Download Psychtoolbox](http://psychtoolbox.org/download)

The code used in the supplemental materials for the Bayesian hierarchical approach employs the following software packages and libraries:
- **Download WinBUGS:** [WinBUGS]( https://www.mrc-bsu.cam.ac.uk/software)
- **WinBUGS/OpenBUGS Manual:** [User Manual](https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/manual14.pdf)
- **R Packages Used:**
  - [R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/index.html)
  - [lme4](https://cran.r-project.org/web/packages/lme4/index.html)

---

# 3. Data Descriptions

### 3.1 Main Data: EPIC Data
The main dataset collected for this project. [Download here]( https://osf.io/jk9nb/?view_only=bfc4c56718334581b267ad2e6f970d74).

### 3.2 Public Data 1: [Robinson and Steyvers' (2023)](https://psycnet.apa.org/manuscript/2023-08265-001.pdf) Flanker Task Data
A public dataset used for analyses requiring larger participant numbers. [Download here](https://osf.io/6hjwv).

### 3.3 Public Data 2: [Hedge et al.'s (2018)](https://link.springer.com/content/pdf/10.3758/s13428-017-0935-1.pdf) Flanker and Stroop Task Data
Another public dataset with larger participant numbers. [Download here](https://osf.io/cwzds).

---

# 4. Script Descriptions
The following scripts are organized into four subfolders based on the datasets and programming tools used: **[EPIC](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/EPIC), [Public1](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/Public1), [Public2](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/Public2),** and **[Bayesian](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/Bayesian)**. The first three subfolders contain Matlab scripts, and the **Bayesian** subfolder includes an R script.

### 4.1 [EPIC Data](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/EPIC)
#### 4.1.1 [Experimental Tasks](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/EPIC/task)
Scripts for running and collecting data from the four inhibitory control tasks used in this study. <br/>
- Flanker Task: `tt_rc_FlankerPracTest.m`
- Prime-Probe Task: `tt_rc_PrimeProbePracTest.m`
- Stroop Task: `tt_rc_StroopTask.m`
- Go/No-Go Task: `tt_rc_GoNoGo.m`

#### 4.1.2 Preprocessing Scripts
Scripts to preprocess EPIC data for plotting violin plots and grand mean plots. <br/>
- `EPIC_preprocess_flanker.m`
- `EPIC_preprocess_primeprobe.m`
- `EPIC_preprocess_Stroop.m`
- `EPIC_violinPlotGrandmeanPlot.m` <br/>

Requires:
- Excel file: `session_numbering.xlsx` (download at [OSF]( https://osf.io/jk9nb/?view_only=bfc4c56718334581b267ad2e6f970d74))
- MATLAB function file: `violinplot.m` (download at [Violinplot-MATLAB](https://github.com/bastibe/Violinplot-Matlab/blob/master/violinplot.m)) <br/>

Outputs: Figure 2, Supplementary Figures 1-4

#### 4.1.3 Rank Order Consistency Across Tasks
A script to rank participants and assess consistency across tasks. <br/>
- `EPIC_rank.m` <br/>

Output: Supplementary Figure 5

#### 4.1.4 Within-Subject Variability Over Time
Scripts to examine within-subject variability over time with extensive repeated testing. <br/>
- `EPIC_practice.m`
- `EPIC_temporal.m` <br/>

Requires:
- `session_numbering.xlsx` <br/>

Outputs: Supplementary Figures 7, 8, 9

#### 4.1.5 Stability Curves
Scripts for generating stability curves. <br/>
- `EPIC_stability_method1.m`
- `EPIC_stability_method2.m`
- `EPIC_stability_Stroop.m`
- `EPIC_stability_ICC.m` <br/>

Requires:
- `session_numbering.xlsx` for all scripts
- MATLAB function file, `ICC.m` ([Download Here](https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)), for `EPIC_stability_ICC.m` <br/>

Outputs: Figure 3, Supplementary Figures 10-15, 17

### 4.2 [Public Data 1](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/Public1)
#### 4.2.1 Stability Curves
A script to generate stability curves using Robinson and Steyvers' (2023) data. <br/>
- `EPIC_Robinson_stability.m` <br/>

Output: Figure 4

#### 4.2.2 Two Model Cases of Small Versus Large Trial Sampling
A script to compare the effects of small and large trial sampling on within- and between-subject variability. <br/>
- `EPIC_Robinson_simulation.m` <br/>

Output: Figure 5

#### 4.2.3 Rank Order Consistency Between the Congruency Effect and Incongruent Trials
A script to rank participants based on their congruency effect and incongruent trial performance to assess consistency. <br/>
- `EPIC_Robinson_rankOrder.m` <br/>

Requires: `ICC.m` <br/>

Outputs: Supplementary Figures 6A, 6D

#### 4.2.4 Drift-Diffusion Modeling
A script to evaluate the reliability of EZ-diffusion model parameters across different trial counts. <br/>
- `EPIC_Robinson_EZdiffusion.m` <br/>

Requires:
- `ezdiffusion.m`
- `ICC.m` <br/>

Outputs: Figure 8, Supplementary Figure 19B

#### 4.2.5 Between-Subject Standard Deviation After Accounting for Trial Noise
A script to correct between-subject standard deviation by adjusting for sampling noise. <br/>
- `EPIC_Robinson_correctingBSsd.m` <br/>

Outputs: Supplementary Figure 18

### 4.3 [Public Data 2](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/Public2)
#### 4.3.1 Simulations to Examine How Within-Subject Variability Contaminates Between-Subject Variability
A script to simulate how varying within-subject variability affects between-subject variability. <br/>
- `EPIC_simulation_variability.m` <br/>

Requires:
- EPIC data
- `session_numbering.xlsx` <br/>

Output: Figure 6

#### 4.3.2 Effect of Trial and Participant Numbers on ICC
A script to plot rank-order consistency as a function of sample size and trial count. <br/>
- `EPIC_simulation_heatmap.m` <br/>

Requires: `ICC.m` <br/>

Output: Figure 7

#### 4.3.3 Rank Order Consistency Between Congruency Effect and Incongruent Trials
Scripts to rank participants based on their congruency effect and incongruent trial performance to assess consistency, with one script incorporating simulations of extended trial sampling. <br/>
- `EPIC_Hedge_rankOrder.m`
- `EPIC_Hedge_simulation_rankOrder.m` <br/>

Requires: `ICC.m` <br/>

Outputs: Supplementary Figures 6B, 6C, 6E, 6F

#### 4.3.4 Drift-Diffusion Model
Scripts to simulate and fit the EZ-diffusion model, examining the effects of sample size and noise on model parameter reliability. <br/>
- `EPIC_Hedge_EZdiffusion.m`
- `EPIC_Hedge_simulation_EZdiffusion_noisesigma.m`
- `EPIC_Hedge_simulation_EZdiffusion_reliability.m` <br/>

Requires:
- `ezdiffusion.m`
- `ICC.m` <br/>

Outputs: Supplementary Figures 19A, 19C, 20

#### 4.3.5 Factor Analysis
Scripts to simulate and conduct confirmatory factor analysis, examining the effects of sample size and noise on factor score reliability. <br/>
- `EPIC_Hedge_simulation_CFA_noisesigma.m`
- `EPIC_Hedge_simulation_CFA_reliability.m` <br/>

Requires: `ICC.m` <br/>

Output: Supplementary Figure 21

### 4.4 [Bayesian Simulation](https://github.com/GrattonLab/LeeSmith_EPIC/tree/main/Bayesian)
This R code simulates 25 replications of reaction time data for a congruency task. Parameters for the simulation were chosen based on estimates from Robinson and Steyvers (2023). The number of trials, subjects, and the ratio of within-subject variance to between-subject variance in the congruency effect can be manipulated at the start of the code.
- `Bayesian_Sim.R` <br/>

1. **Step 1:** A multilevel model using the `lme4` package obtains unbiased estimates of between-subject variability in the congruency effect and trial-level variability in reaction time within subjects. The code generates bias and precision (mean absolute deviation) values for these estimates.

2. **Step 2:** Bayesian estimates of subject-level congruency effects are obtained using WinBUGS. The unbiased estimates from Step 1 are used for variance terms in the priors. The code also generates bias and precision (mean absolute deviation) values for these Bayesian estimates. *Note:* The Bayesian model file ("Cong.txt") should be saved into the working directory for this step.

3. **Step 3:** A more traditional (non-Bayesian) approach estimates subject-level congruency effects. Initial mean reaction times are calculated for each participant for congruent and incongruent trials (no priors are incorporated, so no shrinkage occurs). The difference between these means provides the congruency effect for each subject, along with bias and precision (mean absolute deviation) values.

Output: Supplementary Figure 22

---

# 5. Contact
For any questions or concerns, please email [leehj@illinois.edu].

Last modified on 02/07/2025

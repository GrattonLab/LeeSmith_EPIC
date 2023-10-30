# LeeSmith_EPIC
Last modified on 10/30/2023 <br/>
Any questions or oddities, email Hyejin J. Lee hl22z@fsu.edu <br/>

# Overview
The following is the list of scripts for replicating the EPIC project results (Lee, Smith, Hauenstein, Dworetsky, Kraus, Dorn, Nee, and Gratton, under review). <br/>
You can download the EPIC data at https://osf.io/jk9nb <br/>

# System Requirements
All Matlab scripts are written and tested in Matlab 2021b. <br/>
See this link for Matlab system requirements: https://www.mathworks.com/support/requirements/matlab-system-requirements.html <br/>
and this link to install Matlab: https://www.mathworks.com/help/install/install-products.html <br/>
The installation can take as long as 30 minutes depending on your computer model. <br/>

The code employed in the supplemental materials for comparing traditional and Bayesian estimates uses the following software packages and code libraries. <br/>
Link to download WinBUGS: https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/2018/04/OpenBUGS323setup.zip <br/>
Link to WinBUGS/OpenBUGS manual: https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/manual14.pdf <br/>
R packages used: R2OpenBUGS (https://cran.r-project.org/web/packages/R2OpenBUGS/index.html) <br/>
lme4 (https://cran.r-project.org/web/packages/lme4/index.html)

# Script descriptions
< Preprocessing Matlab scripts > <br/>
Preprocess the EPIC data to plot violin plots and grand mean plots (Fig. 2 and Supp. Fig. 1-5). <br/>
Require an Excel file, "session_numbering.xlsx" (also download at https://osf.io/jk9nb), and Matlab function, violinplot.m (download at https://github.com/bastibe/Violinplot-Matlab/blob/master/violinplot.m). <br/>
1. EPIC_preprocess_flanker <br/>
2. EPIC_preprocess_primeprobe <br/>
3. EPIC_preprocess_Stroop <br/>
4. EPIC_violinPlotGrandmeanPlot <br/>

< Stability curve Matlab scripts > <br/>
Draw stability curves of the EPIC data using methods 1 and 2 (Fig. 3 and Supp. Fig. 7-9). <br/>
Require "session_numbering.xlsx". <br/>
5. EPIC_stability_method1 <br/>
6. EPIC_stability_method2 <br/>

< Public dataset 1 analysis Matlab scripts > <br/>
Analyze Robinson and Steyvers' (2023) flanker task data to draw stability curves (Fig. 4 and Supp. Fig. 10), <br/>
run simulations of two model cases of small versus large trial sampling (Fig. 5), <br/>
examine rank order consistency between congruency effect and incongruent trial performance (Supp. Fig. 6A), <br/>
and conduct drift-diffusion modeling (Fig. 8; require "ezdiffusion.m", which can be downloaded at https://osf.io/jk9nb). <br/>
Require a Matlab function, ICC.m (download at https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc), <br/>
and shadedErrorBar.m (https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar). <br/>
Download Robinson and Steyvers' data at https://osf.io/6hjwv/ <br/>
7. EPIC_Robinson_stability <br/>
8. EPIC_Robinson_simulation <br/>
9. EPIC_Robinson_rankOrder <br/>
10. EPIC_Robinson_EZdiffusion <br/>

< Variability simulation Matlab script > <br/>
Runs a series of simulations to examine how within-subject variability contaimates between-subject variability (Fig. 6) <br/>
No dataset or files required. <br/>
11. EPIC_simulation_variability <br/>

< Public dataset 2 analysis Matlab scripts > <br/>
Analyze and simulate Hedge et al.'s (2018) flanker and Stroop task data to conduct drift-diffusion modeling (Fig. 8 and Supp. Fig. 11), <br/>
factor analysis (Supp. Fig. 12), <br/>
and rank order consistency (Fig. 7 & Supp. Fig. 6B-C). <br/>
Require "ezdiffusion.m" for script#s 13 and 18 and ICC.m. <br/>
Download Hedge et al.'s data at https://osf.io/cwzds/ <br/>
12. EPIC_Hedge_simulation_EZdiffusion_noisesigma <br/>
13. EPIC_Hedge_simulation_EZdiffusion_reliability [Supp. Fig. 11] <br/>
14. EPIC_Hedge_simulation_CFA_noisesigma <br/>
15. EPIC_Hedge_simulation_CFA_reliability [Supp. Fig. 12] <br/>
16. EPIC_Hedge_simulation_rankOrder [Fig 7. & Supp. Fig. 6C] <br/>
17. EPIC_Hedge_rankOrder [Supp. Fig. 6B] <br/>
18. EPIC_Hedge_EZdiffusion [Fig. 8] <br/>

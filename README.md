# LeeSmith_EPIC
Last modified on 09/29/2023
Any questions or oddities, email Hyejin J. Lee hl22z@fsu.edu

The following is the list of scripts for replicating the EPIC project results (Lee, Smith, Hauenstein, Dwortetsky, Kraus, Dorn, Nee, and Gratton, under review)
Download the EPIC data at https://osf.io/jk9nb

< Preprocessing scripts >
Preprocess the EPIC data to plot violin plots and grand mean plots (Fig. 2 and Supp Fig. 1-5).
Require an Excel file, "session_numbering.xlsx" (also download at https://osf.io/jk9nb), and Matlab function, violinplot.m (download at https://github.com/bastibe/Violinplot-Matlab/blob/master/violinplot.m).
1. EPIC_preprocess_flanker
2. EPIC_preprocess_primeprobe
3. EPIC_preprocess_Stroop
4. EPIC_violinPlotGrandmeanPlot

< Stability curve scripts >
Draw stability curves of the EPIC data using methods 1 and 2 (Fig. 3 and Supp. Fig. 7-9).
Also plot a figure showing before and after linear regression.
Require "session_numbering.xlsx".
5. EPIC_stability_method1
6. EPIC_stability_method2

< Robinson and Steyvers' (2023) data analyses scripts >
Analyze Robinson and Steyvers' flanker task data to draw stability curves (Fig. 4 for RT; Supp. Fig. 10 for accuracy), 
run simulations of two model cases (Fig. 5),
examine rank order consistency (Supp. Fig. 6A),
and conduct drift-diffusion modeling (Fig. 8; require "ezdiffusion.m").
Require a Matlab function, ICC.m (download at https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc),
and shadedErrorBar.m (https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar).
(Download Robinson and Steyvers' data at https://osf.io/6hjwv/)
7. EPIC_Robinson_stability
8. EPIC_Robinson_simulation
9. EPIC_Robinson_rankOrder
10. EPIC_Robinson_EZdiffusion

< Variability simulation scripts >
Run a series of simulations to examine how within-subject variability contaimates between-subject variability (Fig. 6)
No dataset or files required.
11. EPIC_simulation_variability

< Hedge et al.'s (2018) data analyses scripts >
Analyze Hedge et al.'s flanker and Stroop task data to conduct drift-diffusion modeling (empirical: Fig. 8; simulated: Supp. Fig. 11),
factor analysis (Supp. Fig. 12), and rank order consistency (Fig. 7 & Supp. Fig. 6B-C)
Require "ezdiffusion.m" for script#s 13 and 18 and ICC.m.
(Download Hedge et al.'s data at https://osf.io/cwzds/)
12. EPIC_Hedge_simulation_EZdiffusion_noisesigma
13. EPIC_Hedge_simulation_EZdiffusion_reliability [Supp. Fig. 11]
14. EPIC_Hedge_simulation_CFA_noisesigma
15. EPIC_Hedge_simulation_CFA_reliability [Supp. Fig. 12]
16. EPIC_Hedge_simulation_rankOrder [Fig 7. & Supp. Fig. 6C]
17. EPIC_Hedge_rankOrder [Supp. Fig. 6B]
18. EPIC_Hedge_EZdiffusion [Fig. 8]

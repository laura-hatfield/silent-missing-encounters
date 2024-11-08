# silent-missing-encounters

Code to implement the analyses reported in Denham and Hatfield "Silent Missingness in Medicaid TAF Data: Can We Fix What We Can’t See?"

# Code to reproduce the simulation results

These files should run out of the box to run the simulation and reproduce the simulation results presented in the manuscript.
**Please note:** the simulation code runs 1000 replications on 10 cores, which may be quite resource-intensive.
Please personalize this for your own computing environment.  

1. `run_silent_missing_simulation.R`   
    Inputs: `silent_missing_simulation_fns.R`  
    Outputs: `silent_missing_sim_results.rds`  
2. `visualize_simulation_results.R`  
    Inputs: `silent_missing_sim_results.rds`  
    Outputs: Figure 2 in the manuscript  


# Code to create the TAF analysis file

For users with access to the TAF data, we include SAS code to produce the analytical files used in the paper.

1. `step_1_ed_visits.txt`  
    Inputs: `taf_other_services_line_2018`  
    Outputs: `ed_vis`
2. `step_2_assign_plan.txt`  
    Inputs: `taf_demog_elig_mngd_care_2018`  
    Outputs: `assigned_plan`  
3. `step_3_classify_ed.txt`  
    Inputs: `ed_vis`  
    Outputs: `ed_vis2`  
4. `step_4_import_comp.txt`  
    Inputs: `composite.xlsx`  
    Outputs: `comp`  
5. `step_5_collapse_ed.txt`  
    Inputs: `ed_vis2`  
    Outputs: `ed_vis2`  
6. `step_6_final_dataset.txt`  
    Inputs: `taf_demo_elig_base_2018`, `fixed_plan_names.xslx`, `assigned_plan`, `ed_vis2`, `taf_apl_base_2018`, `hq`  
    Outputs: `ed_vis_file`  

This includes all the steps to create the sample of plans, enrollees, and to create the plan-level proxy of data completeness. 
Along the way, a few non-TAF external files are used as input, including `composite.xlsx` and `fixed_plan_names.xlsx`.

# Code to create synthetic data and run our approaches on it

Because many users will not have access to the TAF data, we also include code to create a synthetic data set that resembles the TAF data and run the two approaches on it.

1. `create_synthetic_data.R `  
    Inputs:     
    Outputs: `ed_vis_file`  
2. `analyze_synthetic_data.R`  
    Inputs: `ed_vis_file`    
    Outputs:  
3. `visualize_analysis_results.R`  
    Inputs:   
    Outputs: figure analogous to Figure 3 in the manuscript

The code in `analyze_synthetic_data.R` will also run on the data created by the SAS code on the real data (i.e., `ed_vis_file` from step 6 above). 

rule all:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/Fig2.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/model_info_fig2.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/plotdata_fig2.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/descriptive_table.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/tail_density_plot_fig3.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/qq_plot_and_density_plot_fig3.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/quantiles_suppTable3.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/Fig4.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/model_info_fig4.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/plotdata_fig4.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/Fig5.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/model_info_fig5.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/plotdata_fig5.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/prior_history_of_ptd_supp_fig2.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/prior_history_of_ptd_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/prior_history_of_ptd_plotdata.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_supp_fig3.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_plotdata.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/mother_born_preterm_parity_supp_fig4.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/mother_born_preterm_parity_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/mother_born_preterm_parity_plotdata.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/year_of_birth_supp_fig5b.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/year_of_birth_supp_fig5a.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/parity_maternal_age_supp_fig6.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/parity_maternal_age_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/parity_maternal_age_plotdata.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/maternal_characteristics_descriptive_supp_table1.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/simulation_group_bias_supp_table2.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/ultrasound_supp_fig.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/ultrasound_supp_fig_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/plotdata_ultrasound_supp_fig.csv"

## Data cleaning 
# Script filtering
rule filtering_main:
    input:
        "/mnt/hdd/data/swed/ver/mfr_150202_recoded.csv"
    output:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1.csv"
    params:
        "~/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R",
        "~/Parity_gestational_duration_interaction/scripts/functions/2_renumber_parity_to_parityF.R"
    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/filtering.yaml"
    script:
        "scripts/filtering_script_10feb2022.R"


# Script which creates variables based on info from mfr, and other registries which can be connected to the mfr data.
rule create_variables:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1.csv",
        "/mnt/hdd/data/swed/ver/parents.csv",
        "/mnt/hdd/data/swed/ver/education.csv"
    output:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"
    params:
        "~/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"
    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/filtering.yaml"
    script:
        "scripts/create_variables.R"


## Results
# Descriptive table
rule descriptive_table:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv",
        "/mnt/hdd/data/swed/ver/parents.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/descriptive_table.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/descriptive_table.R"


# Regression GD ~ Parity
rule gd_parity:
    input: 
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/Fig2.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/model_info_fig2.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/plotdata_fig2.csv"

    conda: 
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Parity_affects_gd.R"


# Distribution GD
rule variance_gd:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/density_plot_fig3.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/tail_density_plot_fig3.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/qq_plot_fig3.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/quantiles_suppTable3.csv"
        
    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/variance.R"


# Prio history of preterm delivery
rule prio_ptd:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/Fig4.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/model_info_fig4.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/plotdata_fig4.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/clinical_history.R"



# Maternal sister given birth preterm
rule maternal_sister_ptd:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv",
        "/mnt/hdd/data/swed/ver/parents.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/Fig5.eps",
        "/home/karin/Parity_gestational_duration_interaction/results/model_info_fig5.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/plotdata_fig5.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/family_history_maternal_sisters_given_birth_preterm.R"


## Additional Supplementary Results
# Prio history of spontaneous preterm delivery
rule sup_prio_ptd:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/prior_history_of_ptd_supp_fig2.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/prior_history_of_ptd_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/prior_history_of_ptd_plotdata.csv"
    
    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/clinical_history_only_spontaneous.R"


# Family history with no interacting effect
rule supp_fam_history:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv",
        "/mnt/hdd/data/swed/ver/parents.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_supp_fig3.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_plotdata.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/famliy_history_no_parity_interaction_plot.R"


# Three way interaction: mother born preterm*parity*parity of mother
rule three_way_interaction:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/mother_born_preterm_parity_supp_fig4.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/mother_born_preterm_parity_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/mother_born_preterm_parity_plotdata.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/mother_herself_ptd.R"

# Year of delivery interaction
rule YOD:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/year_of_birth_supp_fig5b.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/year_of_birth_supp_fig5a.png"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/year_of_delivery_interaction.R"


# Maternal age
rule supp_maternal_age:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/parity_maternal_age_supp_fig6.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/parity_maternal_age_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/parity_maternal_age_plotdata.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/parity_risk_of_ptd_maternal_age_grouping_mixed_model.R"


# Additional maternal health characteristics, table
rule supp_table:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/maternal_characteristics_descriptive_supp_table1.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/PTB_by_parity_BMI_unwilling_subfertilitu_smoking_diab_preeclamspia_descriptive_table.R"



# Simulation group effect bias
rule simulation_bias:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/simulation_group_bias_supp_table2.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/simulation_group_effect_bias.R"



# Regression GD ~ Parity, only ultrasound. QQ-plot
rule gd_parity_ultrasound:
    input:
        "/mnt/hdd/common/karin/Parity_gestational_duration_interaction/mfr_150202_recored_filtered_p1_all_variables.csv"

    output:
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/ultrasound_supp_fig.png",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/ultrasound_supp_fig_model_info.csv",
        "/home/karin/Parity_gestational_duration_interaction/results/supplementary/plotdata_ultrasound_supp_fig.csv"

    conda:
        "/home/karin/Parity_gestational_duration_interaction/envs/plots.yaml"

    params:
        "/home/karin/Parity_gestational_duration_interaction/scripts/functions/1_cleaning_modules.R"

    script:
        "/home/karin/Parity_gestational_duration_interaction/scripts/plots/Supplementary/ultrasound.R"


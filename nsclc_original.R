#########################################
# Original analyses + step by step exclusions : NSCLC PATIENTS
# OUTCOME : Overall survival
#######################################

##################################
# Renaming and formatting of variables
##################################

#Load data 
df =  read.xlsx(xlsxFile = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09655-y/MediaObjects/41586_2025_9655_MOESM3_ESM.xlsx",
                sheet = "1a Raw Data")
#Renaming variables with special characters in name
df = df %>%
  rename(status_os = `OS_Code_(1=death,.0=censored)`, 
       delay_os = `OS_from_relevant_ICI_start,_months`,
       group = `Group_(1=vaccine)`)

#Set a few parameters  
ylab_txt = "OS (%)"
time_vector = c(12,24,36)
time_for_label = 36
outcome = "OS"

##################################
# Patients' exclusions
##################################

if(same_time_zero){
  df = df %>%
    filter(Previous_Cycles_of_systemic_therapy_prior_to_ICI_start == 0)
}
if(remove_historic_control){
  df = df %>%
    filter(Treatment_Year__ICI >= 2020)
}

#####################################################
# Estimate survival curves
####################################################

#Set order of labels to have the good colors in the plot
df$group = factor(df$group, levels = c(1,0), labels = label_order)

#Estimate weighted Kaplan-Meier survival curves
fit <- survfit(Surv(delay_os, status_os) ~ group,
               data = df)

plot_survival_curves_initial = ggsurvplot(fit,
                                  data = df,
                                  surv.scale = "percent",
                                  xlab = "Months from ICI start",
                                  ylab = ylab_txt,
                                  xlim = c(0,36),
                                  risk.table = TRUE,
                                  legend = "none",
                                  palette = color_for_plots,
                                  legend.labs = label_order,
                                  censor.size = 2,
                                  tables.theme = theme_cleantable(),
                                  ggtheme = theme_light() +    theme(
                                    plot.title = element_text(size = 12)  # ↓ decrease size here
                                  ),
                                  risk.table.col = "strata",
                                  risk.table.title = "Number at risk",
                                  risk.table.fontsize = 2.5,
                                  risk.table.y.text = TRUE,
                                  title = title_txt,
                                  break.time.by = 12)

#Save survival probabilities estimates + 95% CI as dataframe
surv_initial = data.frame(times = summary(fit)$time,
                          surv = summary(fit)$surv,
                          CI_low = summary(fit)$lower,
                          CI_high = summary(fit)$upper,
                          strata = summary(fit)$strata) %>%
  mutate(strata = gsub("group=","", strata)) %>%
  rename(clone = strata)

#Estimate survival at different times (e.g. one, two, three years)
surv_estimates_initial = NULL
for (time in time_vector){
  surv_estimate = surv_initial %>%
    group_by(clone) %>%
    filter(times <= time) %>%
    filter(times == max(times)) %>%
    select(-times) %>%
    mutate(time = time)
  surv_estimates_initial = bind_rows(surv_estimates_initial, surv_estimate)
}

plot_survival_curves = plot_survival_curves_initial$plot
plot_survival_tables = plot_survival_curves_initial$table
surv_point_estimates = surv_estimates_initial

surv_ci_initial = surv_estimates_initial %>%
  group_by(clone, time) %>%
  mutate(label = paste0(round(100*surv,1), "%",
                        ", 95%CI (",
                        round(100*CI_low,1),
                        " - ",
                        round(100*CI_high,1),
                        ")"
  ))

#Create labels for the point estimates + confidence intervals to add as text on the curves
label_ci_vaccine = surv_ci_initial$label[surv_ci_initial$time == time_for_label & surv_ci_initial$clone == "Vaccine"]
label_ci_no_vaccine = surv_ci_initial$label[surv_ci_initial$time == time_for_label & surv_ci_initial$clone == "No vaccine"]
label_tot_initial = paste0(
  "Estimates for 3-year ", outcome, "\nVaccine: ",
  label_ci_vaccine,
  "\nNo vaccine: ",
  label_ci_no_vaccine
)

###########################################################
# Final tuning of the survival curves plot
##########################################################

#Add the label for 3-year OS
plot_survival_curves = plot_survival_curves +
  annotate(geom = "text", x = 1, y = 0.15, label = label_tot_initial, hjust = 0, size = 3)

#Add error bars at 12, 24, 36 months
plot_survival_curves = plot_survival_curves +
  geom_errorbar(data = surv_ci_initial, inherit.aes = FALSE,
               aes(x = time,
               ymin = CI_low, ymax = CI_high,
               color = clone, group = clone),
               position=position_dodge(width=0.5),
               width = 1)

#Trick to align the risk table and the curves vertically
plot_survival_curves_align <- cowplot::align_plots(plot_survival_curves,
                                                   plot_survival_tables + theme(plot.title = element_text(size=10)) +
                                                     theme(legend.position = "none"),
                                                   align = 'v', axis ='l')
plot_survival_curves_nsclc_initial <- cowplot::plot_grid(plot_survival_curves_align[[1]],
                                                 plot_survival_curves_align[[2]],
                                               ncol = 1, rel_heights = c(6,2))


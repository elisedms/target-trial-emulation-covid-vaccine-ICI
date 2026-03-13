##################################
# Libraries and utils
##################################

###############
#Libraries
library(openxlsx)
library(tidyverse)
library(survival)
library(survminer)
library(cowplot)
library(ggimage)

#Util functions
source("util_functions_simulation.R")
syringe_path <- "syringe.png"

###############
#Set seed
set.seed(28)

###############
#Global parameters

label_order = c("Vaccine", "No vaccine")
#Blue for vaccine - red for no vaccine
color_for_plots = c("#00649A", "#c1121f")

#Survival probabilities for data simulation
times <- c(0, 25, 50, 75, 100, 175, 250)
surv  <- c(1, 0.8, 0.70, 0.65, 0.60, 0.50, 0.40)

##For vaccination
#Probability of treatment (at the end of the grace period)
p_treatment = 0.5
#Hazard ratio for treatment versus control
hr_vaccination_susceptibles = 3
hr_vaccination_not_susceptibles = 0.8
#Length of the grace period
length_grace_period = 100

#For susceptibles versus not-susceptibles
p_susceptibles = 0.5
hr_susceptibles = 2

#Number of patients in simulation
N = 10000

###############
#Data simulation

#First, we simulate vaccination times
times_vaccination = times[times <= length_grace_period]
surv_vaccination = seq(from = 1, to = p_treatment, length.out = length(times_vaccination))
dat_vaccination <- simulate_piecewise_treatment(
  times = times_vaccination,
  surv  = surv_vaccination,
  N     = N,
  tau   = length_grace_period
)
dat_vaccination = dat_vaccination %>%
  rename(time_vaccination = time,
         vaccine_status = event) %>%
  #To disentangle the effect of immortal-time bias
  #We assume that vaccinated patients are all vaccinated at time 0
  mutate(time_vaccination = ifelse(vaccine_status == 1,0,Inf))

#Second, we add an indicator for frailty : susceptibles or not susceptibles?
dat_vaccination$susceptibility_status = rbinom(N,size = 1, prob = p_susceptibles)

# Third, we simulate survival times
# (based on hazard ratio for treatment and for susceptibles versus non-susceptibles)
dat_sim_tot <- simulate_piecewise_survival(
  dat_vaccination, hr_vaccination_susceptibles, hr_vaccination_not_susceptibles,
  hr_susceptibles,
  times, surv,
  N,
  tau = max(times)
)

###################################
# Left plot : hypothetical study

# Whole population
dat_sim_tot$strata = case_when(
  dat_sim_tot$vaccine_status == 0 & dat_sim_tot$susceptibility_status == 0 ~ "Unvaccinated - not susceptible",
  dat_sim_tot$vaccine_status == 0 & dat_sim_tot$susceptibility_status == 1 ~ "Unvaccinated - susceptible",
  dat_sim_tot$vaccine_status == 1 & dat_sim_tot$susceptibility_status == 0 ~ "Vaccinated - not susceptible",
  dat_sim_tot$vaccine_status == 1 & dat_sim_tot$susceptibility_status == 1 ~ "Vaccinated - susceptible"
)
label_order_susceptible_group = c("Vaccinated - not susceptible",
                "Vaccinated - susceptible",
                "Unvaccinated - not susceptible",
                "Unvaccinated - susceptible"

)
dat_sim_tot$strata = factor(dat_sim_tot$strata, levels = label_order_susceptible_group)

fit_wp <- survfit(Surv(time, event) ~ strata, data = dat_sim_tot)
survival_curve_wp = ggsurvplot(fit_wp,
                               data = dat_sim_tot,
                               surv.scale = "percent",
                               xlab = "Days from baseline",
                               ylab = "Overall survival",
                               xlim = c(0,250),
                               risk.table = FALSE,
                               legend = "none",
                               palette = rep(color_for_plots,each = 2), 
                               linetype = rep(c("dotted", "twodash"), times = 2),
                               legend.labs = label_order_susceptible_group,
                               legend.title = "",
                               conf.int = FALSE,
                               censor.size = 0,
                               tables.theme = theme_cleantable(),
                               ggtheme = theme_minimal())

plot_setting <- survival_curve_wp$plot + 
  ggtitle("Survival-induced selection bias",
          subtitle = "Hypothetical study")

##############################
# Middle plot : survival curves whole population

#Plot sample patient is replaced by the survival curves per treatment group on average
data_sim_tot_middle_plot = dat_sim_tot %>%
  mutate(vaccine_status = ifelse(vaccine_status == 1, "Vaccine", "No vaccine")) %>%
  mutate(vaccine_status = factor(vaccine_status, levels = label_order))
fit_wp <- survfit(Surv(time, event) ~ vaccine_status, data = data_sim_tot_middle_plot)
plot_sample_patient = ggsurvplot(fit_wp,
                               data = data_sim_tot_middle_plot,
                               surv.scale = "percent",
                               xlab = "Days from baseline",
                               ylab = "Overall survival",
                               xlim = c(0,250),
                               risk.table = FALSE,
                               legend = "none",
                               palette = color_for_plots, 
                               legend.title = "",
                               conf.int = FALSE,
                               censor.size = 0,
                               tables.theme = theme_cleantable(),
                               ggtheme = theme_minimal())

plot_sample_patient = plot_sample_patient$plot + 
  ggtitle("",
          subtitle = "Survival in the whole population")

##############################
# Right plot : survival curves per treatment group

# By treatment group (with the colors)
dat_sim_alternative_time_zero = dat_sim_tot %>%
  #We do not count events occurring after 100 days  
  mutate(event = ifelse(event == 1 & time <= 100, 0, event)) %>%
  mutate(vaccine_status = ifelse(vaccine_status == 1, "Vaccine", "No vaccine")) %>%
  mutate(vaccine_status = factor(vaccine_status, levels = label_order))

fit_per_treatment_group <- survfit(Surv(time, event) ~ vaccine_status, data = dat_sim_alternative_time_zero)
survival_curve_per_treatment_group = ggsurvplot(fit_per_treatment_group,
                                                data = dat_sim_alternative_time_zero,
                                                surv.scale = "percent",
                                                xlab = "Days from baseline",
                                                ylab = "Overall survival",
                                                xlim = c(0,250),
                                                risk.table = FALSE,
                                                legend = "none",
                                                palette = color_for_plots, 
                                                legend.labs = c("Vaccine","No vaccine"),
                                                conf.int = FALSE,
                                                censor.size = 0,
                                                tables.theme = theme_cleantable(),
                                                ggtheme = theme_minimal())


###################################
#Create legend

# Order
order_levels <- c(
  "Vaccine group",
  "Control group"
)

# Dummy data
df <- data.frame(
  x = 1,
  y = 1,
  group = factor(order_levels, levels = order_levels)
)

plot_legend = ggplot(df, aes(x, y)) +
  geom_segment(
    data = df[df$group %in% c("Vaccine group", "Control group"), ],
    aes(x = 0, xend = 1, yend = y,
        color = group,),
    linewidth = 1
  ) +
  scale_color_manual(
    breaks = order_levels,
    values = c(
      "Vaccine group" =  color_for_plots[1],
      "Control group" =  color_for_plots[2]
    ),
    na.translate = FALSE
  ) +
  guides(
    color = guide_legend(nrow = 1)#,
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
legend = cowplot::get_legend(plot_legend)

###################################
#Plot assembly
plot_tot_rawC = cowplot::plot_grid(plot_setting,
                                   plot_sample_patient,
                                   survival_curve_per_treatment_group$plot + ggtitle("", 
                                                                                     subtitle = "Removal of events occuring before 100 days"),
                                   nrow = 1)
plot_tot_rawC = cowplot::plot_grid(plot_tot_rawC, legend, ncol = 1, rel_heights = c(10,1))
plot_tot_rawC

#Save plot
pdf("Figure2_rawC.pdf",
    width = unit(1.2*10,"cm"),
    height = unit(1.1*1.2*3, "cm"))
  plot_tot_rawC
dev.off()

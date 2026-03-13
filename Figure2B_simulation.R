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

###############
#Set seed
set.seed(28)

#Util functions
source("util_functions_simulation.R")
syringe_path <- "syringe.png" 

###############
#Global parameters

label_order = c("Vaccine", "No vaccine")
#Blue for vaccine - red for no vaccine
color_for_plots = c("#00649A", "#c1121f")

#Survival probabilities for data simulation
times_surv <- c(0, 25, 50, 75, 100, 175, 250, 350)
prob_surv  <- c(1, 0.95, 0.90, 0.85, 0.80, 0.65, 0.50, 0.35)

##For vaccination
#Probability of treatment
p_treatment = 0.5 #Treatment is protective
#Hazard ratio for treatment versus control
hr_vaccination_susceptibles = 2
hr_vaccination_not_susceptibles = 1.15

#Length of the grace period
length_grace_period = 100

#For susceptibles versus not-susceptibles
p_susceptibles = 0.5
#hr_susceptibles = 10
hr_susceptibles = 5

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
  #For now we imagine they are all vaccinated at -100 days
  mutate(time_vaccination = ifelse(vaccine_status == 0, Inf, 0))

#Second, we add an indicator for frailty : susceptibles or not susceptibles?
dat_vaccination$susceptibility_status = rbinom(N,size = 1, prob = p_susceptibles)

# Third, we simulate survival times
# (based on hazard ratio for treatment and for susceptibles versus non-susceptibles)
dat_sim_tot <- simulate_piecewise_survival(
  dat_vaccination,
  hr_vaccination_susceptibles, 
  hr_vaccination_not_susceptibles, 
  hr_susceptibles,
  times = times_surv, 
  surv = prob_surv,
  N,
  tau = max(times_surv)
)


###################################
# Left plot : hypothetical study

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
                                  xlim = c(0,350),
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
  scale_x_continuous(
    labels = function(x) x - 100
  ) +
  ggtitle("Prevalent-user bias",
          subtitle = "Hypothetical study") +
  geom_vline(xintercept = 100, color = "#868686", linetype = "dotted") +
  annotate(geom = "text", x = 105, y = 1, hjust = 0, label = "Eligibility")

###################################
# Middle plot : sample of 10 patients

# Sample of 10 individuals : 5 with treatment, 5 without treatment
# Draw a line plot for these individuals

dat_sim_tot$vaccine_status = ifelse(dat_sim_tot$vaccine_status == 1,
                                           "Vaccine", "No vaccine")
dat_sim_tot$vaccine_status = factor(dat_sim_tot$vaccine_status, levels = label_order)

dat_sim_tot_sample = dat_sim_tot %>% 
  group_by(vaccine_status, event) %>% 
  slice_sample(n = 3) %>% 
  group_by(vaccine_status) %>%
  slice_sample(n = 5) %>% 
  ungroup() %>%
  arrange(desc(vaccine_status)) %>%
  mutate(sample_number = row_number())

dat_plot <- dat_sim_tot_sample %>%
  mutate(syringe = if_else(vaccine_status == "Vaccine", syringe_path, NA_character_)) %>%
  mutate(time_vaccination = ifelse(vaccine_status == "No vaccine", 0, time_vaccination)) %>%
  mutate(vaccine_eligibility_status = case_when(
    time < 100 ~ "uneligible",
    TRUE ~ vaccine_status
  )) %>%
  mutate(vaccine_eligibility_status = factor(vaccine_eligibility_status,
                                             levels = c("uneligible", "Vaccine","No vaccine")))
offset_y <- 0.42
offset_x <- -3

plot_sample_patient = ggplot(dat_plot,
       aes(y = sample_number, x = 0, xend = time, yend = sample_number, color = vaccine_eligibility_status)) +
  geom_segment(linewidth = 1, linetype = "dotted", aes(y = sample_number, x = 0, xend = time_vaccination, yend = sample_number)) +
  geom_segment(linewidth = 1,aes(y = sample_number, x = time_vaccination, xend = time, yend = sample_number)) +
  # Add syringe icons at vaccination time for vaccinated rows
  geom_image(
    data = ~ dplyr::filter(.x, vaccine_status == "Vaccine"),
    aes(x = time_vaccination + offset_x, y = sample_number + offset_y, 
        image = syringe),
    size = 0.1, 
    asp = 1
  ) +
  geom_point(
    data = ~ dplyr::filter(.x, event == 1),
    aes(x = time, y = sample_number),
    shape = 23,        
    size  = 1,
    fill  = "white",  
    stroke = 1.2   
  ) +
  scale_color_manual(values = c("#868686", color_for_plots)) +
  theme_minimal() +
  labs(x = "Days from baseline", y = "Patient ID") +
  theme(legend.position = "none") +
  scale_x_continuous(
    labels = function(x) x - 100
  ) +
  scale_y_continuous(
    breaks = 1:10,
    limits = c(0.5,10.9)
  ) +
  geom_vline(xintercept = 100, color = "#868686", linetype = "dotted") +
  annotate(geom = "text", x = 105, y = 10.8, hjust = 0, label = "Eligibility")

plot_sample_patient

##############################
# Right plot : survival curves per treatment group

# By treatment group (with the colors)

dat_sim_eligibile = dat_sim_tot %>%
  filter(time > 100) %>%
  mutate(time = time - 100,
         time_vaccination = time_vaccination - 100)
dim(dat_sim_eligibile)
summary(factor(dat_sim_eligibile$vaccine_status))
summary(factor(dat_sim_eligibile$susceptibility_status))

fit_per_treatment_group <- survfit(Surv(time, event) ~ vaccine_status, data = dat_sim_eligibile)
survival_curve_per_treatment_group = ggsurvplot(fit_per_treatment_group,
                               data = dat_sim_eligibile,
                               surv.scale = "percent",
                               xlab = "Days from baseline",
                               ylab = "Overall survival",
                               xlim = c(0,250),
                               risk.table = FALSE,
                               legend = "none",
                               palette = color_for_plots, 
                               conf.int = FALSE,
                               censor.size = 0,
                               tables.theme = theme_cleantable(),
                               ggtheme = theme_minimal())

print(survival_curve_per_treatment_group)

###################################
#Create legend

# Order
order_levels <- c(
  "Vaccine group",
  "Control group",
  "Excluded patient",
  "Outcome",
  "Vaccination date"
)

# Dummy data
df <- data.frame(
  x = 1,
  y = 1,
  group = factor(order_levels, levels = order_levels)
)

# Custom legend key to draw the syringe
draw_key_syringe <- function(data, params, size) {
  grid::rasterGrob(
    png::readPNG(syringe_path),
    width = unit(1, "npc"),
    height = unit(1, "npc")
  )
}

plot_legend = ggplot(df, aes(x, y)) +
  geom_segment(
    data = df[df$group %in% c("Vaccine group", "Control group", "Excluded patient"), ],
    aes(x = 0, xend = 1, yend = y,
        color = group,),
    linewidth = 1
  ) +
  geom_point(
    data = df[df$group == "Outcome", ],
    aes(shape = group),
    size = 1,
    fill = "white",
    color = "#868686"
  ) +
  geom_point(
    data = df[df$group == "Vaccination date", ],
    aes(shape = group),
    size = 3,
    key_glyph = draw_key_syringe
  ) +
  scale_color_manual(
    breaks = order_levels,
    values = c(
      "Vaccine group" =  color_for_plots[1],
      "Control group" =  color_for_plots[2],
      "Excluded patient" =  "#868686"
    ),
    na.translate = FALSE
  ) +
  scale_shape_manual(
    breaks = order_levels,
    values = c(
      "Outcome" = 23,
      "Vaccination date" = NA
    ),
    na.translate = FALSE
  ) +
  guides(
    color = guide_legend(nrow = 1),
    #linetype = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
legend = cowplot::get_legend(plot_legend)

###################################
#Plot assembly
plot_tot_rawB = cowplot::plot_grid(plot_setting,
                                   plot_sample_patient + 
                                         ggtitle("", 
                                               subtitle = "Sample of 10 patients"),
                                     
                                   survival_curve_per_treatment_group$plot + ggtitle("", 
                                                                                     subtitle = "Survival in eligible patients by vaccine group"),
                                   nrow = 1)
plot_tot_rawB = cowplot::plot_grid(plot_tot_rawB, legend, ncol = 1, rel_heights = c(10,1))
plot_tot_rawB

#Save plot
pdf("Figure2_rawB.pdf",
    width = unit(1.2*10,"cm"),
    height = unit(1.1*1.2*3, "cm"))
    plot_tot_rawB
dev.off()

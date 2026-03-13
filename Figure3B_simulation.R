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
injection_path = "medical.png"

###############
#Set seed
set.seed(28)

###############
#Global parameters

label_order = c("Vaccine", "No vaccine")
#Blue for vaccine - red for no vaccine
color_for_plots = c("#00649A", "#c1121f")

#Survival probabilities for data simulation
times <- c(0, 48*25/250, 48*50/250, 48*75/250, 48*100/250, 48*175/250, 50*250/250)
surv  <- c(1, 0.8, 0.70, 0.65, 0.60, 0.50, 0.40)

##For vaccination
#Probability of treatment (at the end of the grace period)
p_treatment = 0.5
#Hazard ratio for treatment versus control
hr_vaccination_susceptibles = 1
hr_vaccination_not_susceptibles = 1
#Length of the grace period
length_grace_period = 48

#For susceptibles versus not-susceptibles
p_susceptibles = 0.5
hr_susceptibles = 1

#Number of patients in simulation
N = 10000

###############
#Data simulation

#First, we simulate vaccination times
possible_dates_vaccination = seq(0,50,9)

dat_vaccination <- simulate_piecewise_treatment(
  times = times_vaccination,
  surv  = surv_vaccination,
  N     = N,
  tau   = length_grace_period
)

dat_vaccination <- data.frame(
  id    = seq_len(N),
  vaccine_status = rbinom(N,size = 1, prob = p_treatment)
)

dat_vaccination <- dat_vaccination %>%
  rowwise() %>%
  mutate(time_vaccination = ifelse(vaccine_status == 0, 48, sample(possible_dates_vaccination, 1)))

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

# Next, change baseline for vaccinated patients
dat_sim_tot <- dat_sim_tot %>%
  mutate(time_zero = case_when(
    vaccine_status == 0 ~ 0, #Date of ICI initiation
    vaccine_status == 1 ~ 9 * round(time_vaccination / 9), #Closest ICI treatment line
  )) %>%
  #Correction : if patients die between vaccination and last multiple of 9 months, we assign to the injection before
  mutate(time_zero = ifelse(time_zero > time,  time_zero - 9, time_zero)) %>%
  mutate(time_updated = time - time_zero)

###################################
# Left plot : hypothetical study

# Whole population
fit_wp <- survfit(Surv(time, event) ~ 1, data = dat_sim_tot)
survival_curve_wp = ggsurvplot(fit_wp,
                              data = dat_sim_tot,
                              surv.scale = "percent",
                              xlab = "Months from first ICI",
                              ylab = "Overall survival",
                              xlim = c(0,48),
                              ylim = c(0,1.05),
                              risk.table = FALSE,
                              legend = "none",
                              palette = "black",
                              conf.int = FALSE,
                              censor.size = 0,
                              tables.theme = theme_cleantable(),
                              ggtheme = theme_minimal())

p <- survival_curve_wp$plot

plot_setting <- p +
  ggtitle("Bias due to differential time zero",
          subtitle = "Hypothetical study") +
  
  geom_point(data = data.frame(x = seq(from = 0, to = 48, by = 9) ,
                                 y = 0),
               inherit.aes = FALSE,
               aes(x =x , y = y),
               color = color_for_plots[1]) +
  geom_segment(data = data.frame(x = seq(from = 0, to = 48, by = 9),
                                 xend = seq(from = 0, to = 48, by = 9),
                                 y = 0, yend = 1),
               inherit.aes = FALSE,
               aes(x =x , xend = xend, y = y, yend = yend),
               color = "#868686", linetype = "dotted") +
  geom_image(
    data = ~ data.frame(x = seq(from = 0, to = 48, by = 9), y = 1.05, injection = injection_path),
    aes(x = x, y = y, 
        image = injection),
    size = 0.05,
    asp = 1,
    inherit.aes = FALSE
  )

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
  mutate(time_vaccination = ifelse(vaccine_status == "No vaccine", 0, time_vaccination))
offset_y <- 0.42
offset_x <- -3*48/250

plot_sample_patient = ggplot(dat_plot,
       aes(y = sample_number, x = 0, xend = time, yend = sample_number, color = vaccine_status)) +
  geom_segment(data = data.frame(x = seq(from = 0, to = 48, by = 9),
                                 xend = seq(from = 0, to = 48, by = 9),
                                 y = 0.9, yend = 11),
                 inherit.aes = FALSE, 
                aes(x =x , xend = xend, y = y, yend = yend),
                color = "#868686", linetype = "dotted") +
  geom_image(
    data = ~ data.frame(x = seq(from = 0, to = 48, by = 9), y = 11.3, injection = injection_path),
    aes(x = x, y = y, 
        image = injection),
    size = 0.05,
    asp = 1,
    inherit.aes = FALSE
  ) +
  geom_segment(linewidth = 1,aes(y = sample_number, x = time_zero, xend = time, yend = sample_number)) +
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
  scale_color_manual(values = color_for_plots) +
  theme_minimal() +
  labs(x = "Months from first ICI", y = "Patient ID") +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = 1:10,
    limits = c(0.5,11.3)
  )

##############################
# Right plot : survival curves per treatment group

# By treatment group (with the colors)
dat_sim_tot$vaccine_status = factor(dat_sim_tot$vaccine_status, levels = label_order)
fit_per_treatment_group <- survfit(Surv(time_updated, event) ~ vaccine_status, data = dat_sim_tot)
survival_curve_per_treatment_group = ggsurvplot(fit_per_treatment_group,
                               data = dat_sim_tot,
                               surv.scale = "percent",
                               xlab = "Months from baseline",
                               ylab = "Overall survival",
                               xlim = c(0,48),
                               risk.table = FALSE,
                               legend = "none",
                               palette = color_for_plots, 
                               #legend.labs = label_order,
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
    data = df[df$group %in% c("Vaccine group", "Control group"), ],
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
      "Control group" =  color_for_plots[2]
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
plot_tot_rawE = cowplot::plot_grid(plot_setting,
                                   plot_sample_patient + 
                                         ggtitle("", 
                                                 subtitle = "Sample of 10 patients"),
                                   survival_curve_per_treatment_group$plot + ggtitle("", 
                                   subtitle = "Survival in eligible patients by vaccine group"),
                                   nrow = 1)
plot_tot_rawE = cowplot::plot_grid(plot_tot_rawE, legend, ncol = 1, rel_heights = c(10,1))
plot_tot_rawE

#Save plot
pdf("Figure2_rawE.pdf",
    width = unit(1.2*10,"cm"),
    height = unit(1.1*1.2*3, "cm"))
    plot_tot_rawE
dev.off()

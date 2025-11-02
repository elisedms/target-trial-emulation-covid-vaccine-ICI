##########################################################################
# This script is a util file : it provides function used in the main codes
##########################################################################

##########################################################################
# Clone and censor
##########################################################################
clone_censor <- function(df_restricted,
                         end_grace_period_month = 100/30.436875){
  
  ########
  #List of parameters :
  #df is the dataframe
  #end_grace_period_month is the length of the grace period (same unit than delay_os - here in months) - one month is 30.436875 days on average
  ########
  
  df_restricted_clone = 
    #Cloning 
    bind_rows(df_restricted %>% mutate(clone = "Vaccine"),
              df_restricted %>% mutate(clone = "No vaccine")) %>%
    #Artificial censoring
    mutate(status_os = case_when(
      #Clones in the vaccine group who did not receive the vaccine at the end of the grace period are censored at the end of the grace period
      clone == "Vaccine" & group == "Observed no vaccine" & delay_os > end_grace_period_month  ~ 0,
      #Clones in the no vaccine group who did receive the vaccine should be censored at the time they receive the vaccine
      #However, as we do not have the precise date of vaccination, we censored them at the end of the grace period, as discussed in the target trial emulation protocol
      clone == "No vaccine" & group == "Observed vaccine" & delay_os > end_grace_period_month  ~ 0,
      #Otherwise clones are not artificially censored
      TRUE ~ status_os
    )) %>%
    mutate(delay_os = case_when(
      clone == "Vaccine" & group == "Observed no vaccine" & delay_os > end_grace_period_month  ~ end_grace_period_month,
      clone == "No vaccine" & group == "Observed vaccine" & delay_os > end_grace_period_month ~ end_grace_period_month,
      TRUE ~ delay_os
    )) %>%
    #We create a boolean variable indicator of if the clone is censored
    #(useful to restrict the data to observed treatment patterns -- only for plots and descriptive counts)
    mutate(clone_boolean = case_when(
      clone == "Vaccine" & group == "Observed no vaccine" ~ TRUE,
      clone == "No vaccine" & group == "Observed vaccine" ~ TRUE,
      TRUE ~ FALSE
    ))
  
  return(df_restricted_clone)
}

##########################################################################
# Estimate inverse probability of treatment weights
##########################################################################

weight <- function(df_restricted_clone,
                   variables_used_for_adjustment,
                   print = TRUE,
                   color_for_plots = NULL,
                   end_grace_period_month = 100/30.436875){
  
  ########
  #List of parameters :
  #df_restricted_clone is the dataframe (after cloning)
  #variables_used_for_adjustment : list of variables used for adjustment
  #print : TRUE if the quality check plots for propensity score (PS) overlap and weights distribution should be plotted. Typically set to FALSE during bootstrap
  #color_for_plots : vector of colors used for the plot
  #end_grace_period_month is the length of the grace period (same unit than delay_os - here in months) - one month is 30.436875 days on average
  ########
  
  ################################################
  #Fit logistic model to estimate propensity scores
  
  #Response variable should be numerical
  df_restricted_clone$group_num = ifelse(df_restricted_clone$group == "Observed vaccine",1,0) 
  formula_variable = paste0(variables_used_for_adjustment, collapse = "+")
  #The propensity score model is fitted on the uncloned data
  ps_model = glm(formula = paste0("group_num ~ ", formula_variable),
                 data = df_restricted_clone %>% filter(!clone_boolean),
                 family = quasibinomial())
  
  if(print){
    #Print the model summary
    print("Fitted propensity score model")
    print(summary(ps_model))
  }
  
  ################################################
  #Estimate the propensity score on the cloned dataset
  
  df_restricted_clone$ps = predict(ps_model,
                                   newdata = df_restricted_clone,
                                   type = "response")
  
  #Relevel variable group to have the factor in the right order for the plot
  df_restricted_clone$group = factor(df_restricted_clone$group,
                                     levels = c("Observed vaccine", "Observed no vaccine"),
                                     labels = label_order)
  df_restricted_no_clone = df_restricted_clone %>% filter(!clone_boolean)
  
  #Check propensity score overlap
  if(print){
    plot_overlap = ggplot(data = df_restricted_no_clone,
           aes(x=ps, group = group, fill = group)) +
      geom_density(alpha = 0.8) +
      scale_fill_manual(values = color_for_plots) + 
      theme_bw() +
      theme(legend.position = "bottom") + 
      labs(x = "Propensity score", y = "Density", fill = "Observed treatment group")
  }
  
  #Estimate the weights based on the propensity score
  df_restricted_clone$weight = case_when(
    df_restricted_clone$group_num == 1 ~ 1/df_restricted_clone$ps,
    df_restricted_clone$group_num == 0 ~ 1/(1-df_restricted_clone$ps)
  )
  
  #Check distribution of weights
  if(print){
    plot_weight_distribution =   ggplot(data = df_restricted_clone %>% filter(!clone_boolean),
                                        aes(group = group, x = group, fill = group, y = weight)) + 
      geom_boxplot() +
      scale_fill_manual(values = color_for_plots) + 
      theme_bw() +
      theme(legend.position = "none") + 
      labs(x = "Observed treatment", y = "Weight distribution")
  }
  
  #Before the end of the grace period -- nobody can be censored so weighs are 1 for everyone
  #At the end of the grace period uncensored patients are up-weighted with the weight computed above
  df_restricted_clone_long =
    bind_rows(df_restricted_clone %>% mutate(tstart = 0, tstop = end_grace_period_month, ipcw = 1),
              df_restricted_clone %>% mutate(tstart = end_grace_period_month, tstop = delay_os, ipcw = weight)) %>%
    #For patients who experience the event during the grace period
    filter(tstart < delay_os) %>%
    mutate(tstop = pmin(tstop, delay_os)) %>%
    #Event only at the final row
    mutate(status_os = ifelse(tstop == delay_os, status_os, 0))
  
  if(print){
    return(list(df_restricted_clone_long, df_restricted_clone, plot_overlap, plot_weight_distribution))
  }else{
    return(list(df_restricted_clone_long, df_restricted_clone))
  }
  
}

##########################################################################
# Estimate weighted Kaplan-Meier survival curves and 3-year OS
##########################################################################

estimate_survival <- function(df_restricted_clone_long,
                              df_restricted_clone,
                   color_for_plots = NULL,
                   ylab_txt = "OS (%)",
                   print = TRUE,
                   time_vector = c(36),
                   title_txt = "NSCLC + ICI"
                   ){
  
  ########
  #List of parameters :
  #df_restricted_long_clone is the dataframe in long format after cloning (one row during grace period -- one row after grace period)
  #df_restricted_clone is the dataframe (after cloning)
  #color_for_plots : vector of colors used for the weights
  #ylab_txt : lab for the y axis of the survival curves
  #print : TRUE if the survival plots should be plotted (set to FALSE for bootstrapping)
  #time_vecotr : Vector of times when estimates of survival are given (same unit that delay_os)
  #title_txt : title for the plot
  ########
  
  #Estimate weighted Kaplan-Meier survival curves on the weighted population
  fit <- survfit(Surv(tstart, tstop, status_os) ~ clone,
                 data = df_restricted_clone_long,
                 weights = df_restricted_clone_long$ipcw)
  
  #Plot the survival curves
  if(print){
    
    plot_survival_curves = ggsurvplot(fit,
                                      data = df_restricted_clone_long,
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
                                      ggtheme = theme_light(),
                                      title = title_txt)
    
    #For the risk tables we use cloned but unweighted estimates
    fit_unweighted <- survfit(Surv(delay_os, status_os) ~ clone,
                   data = df_restricted_clone)
    
    plot_survival_curves_raw = ggsurvplot(fit_unweighted,
                                      data = df_restricted_clone,
                                      surv.scale = "percent",
                                      xlim = c(0,36),
                                      risk.table = TRUE,
                                      legend = "none",
                                      palette = color_for_plots,
                                      risk.table.col = "strata",
                                      risk.table.title = "Number at risk (unweighted)",
                                      risk.table.fontsize = 2.5,
                                      legend.labs = label_order,
                                      risk.table.y.text = TRUE,
                                      tables.theme = theme_cleantable(),
                                      ggtheme = theme_light(),
                                      title = title_txt)
    
    
  }
  
  #Save survival probabilities as dataframe
  surv_weighted = data.frame(times = summary(fit)$time,
                                     surv = summary(fit)$surv,
                                     strata = summary(fit)$strata) %>%
    mutate(strata = gsub("clone=","", strata)) %>%
    rename(clone = strata)
  
  #Estimate survival at different times (e.g. one, two, three years)
  surv_estimates = NULL
  for (time in time_vector){
    surv_estimate = surv_weighted %>%
      group_by(clone) %>%
      filter(times <= time) %>%
      filter(times == max(times)) %>%
      select(-times) %>%
      mutate(time = time)
    surv_estimates = bind_rows(surv_estimates, surv_estimate)
  }
  
  if(print){
    return(list(plot_survival_curves$plot,plot_survival_curves_raw$table,
                surv_estimates))
  }else{
    return(list(surv_estimates))
  }
  
}

##########################################################################
# Estimate 95% confidence intervals from bootstrap samples
##########################################################################

estimate_ci <- function(surv_point_estimates, surv_estimates_bootstrap,
                        outcome = "OS", time_for_label = 36){

  ########
  #List of parameters :
  #surv_point_estimates : data frame with point estimates of survival on original data
  #surv_estimates_bootstrap : data frame with point estimates of survival on bootstrapped data
  #outcome is a string for the outcome (default :  "OS") -- used for label_tot (text for the survival plot)
  ########
  
  #
  surv_ci= surv_estimates_bootstrap %>%
    group_by(clone, time) %>%
    summarise(CI_low = unname(quantile(surv,probs = c(0.025)))[1],
              CI_high = unname(quantile(surv,probs = c(0.975)))[1]) %>%
    left_join(surv_point_estimates) %>%
    mutate(label = paste0(round(100*surv,1), "%",
           ", 95%CI (",
           round(100*CI_low,1),
           " - ",
           round(100*CI_high,1),
           ")"
     ))
  
  #Create labels for the point estimates + confidence intervals to add as text on the curves
  label_ci_vaccine = surv_ci$label[surv_ci$time == time_for_label & surv_ci$clone == "Vaccine"]
  label_ci_no_vaccine = surv_ci$label[surv_ci$time == time_for_label & surv_ci$clone == "No vaccine"]
  label_tot = paste0(
    "Estimates for 3-year ", outcome, "\nVaccine: ",
    label_ci_vaccine,
    "\nNo vaccine: ",
    label_ci_no_vaccine
  )
  
  #Difference between arms at time_for_label
  surv_estimate_diff = surv_estimates_bootstrap$surv[surv_estimates_bootstrap$time == time_for_label & surv_estimates_bootstrap$clone == "Vaccine"] - 
                       surv_estimates_bootstrap$surv[surv_estimates_bootstrap$time == time_for_label & surv_estimates_bootstrap$clone == "No vaccine"] 
  surv_estimate_diff_low = unname(quantile(surv_estimate_diff,probs = c(0.025)))[1]
  surv_estimate_diff_high = unname(quantile(surv_estimate_diff,probs = c(0.975)))[1]
  surv_point_estimate_vaccine_diff = surv_point_estimates$surv[surv_point_estimates$time == time_for_label & surv_point_estimates$clone == "Vaccine"]-
                                     surv_point_estimates$surv[surv_point_estimates$time == time_for_label & surv_point_estimates$clone == "No vaccine"]
  
  label_ci_diff = paste0(
    round(100*surv_point_estimate_vaccine_diff,1), "%",
    ", 95%CI (",
    round(100*surv_estimate_diff_low,1),
    " - ",
    round(100*surv_estimate_diff_high,1),
    ")"
  )
  
  return(list(
    surv_ci,
    label_ci_diff,
    label_tot
  ))
  
}

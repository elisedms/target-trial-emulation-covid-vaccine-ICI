##############################################################
#Util functions for data simulation
##############################################################

##############################################################
#Simulate survival status - with piecewise constant hazard function
##############################################################
simulate_piecewise_survival <- function(dat_vaccination, 
                                        hr_vaccination_susceptibles, 
                                        hr_vaccination_not_susceptibles, 
                                        hr_susceptibles,
                                        times, surv,
                                        N,
                                        tau = max(times)) {
  
  
  # Interval lengths and hazards
  dt <- diff(times)
  lambda <- log(surv[-length(surv)] / surv[-1]) / dt
  
  #For each unit, depending on vaccination status and time
  #the hazard gets multiplied by hr_vaccination after initiation of treatment
  #And by the hr_susceptibles for susceptible patients
  df_hazard = data.frame(time = times %>% head(-1), dt = dt, lambda = lambda)
  data_sim = dat_vaccination %>%
    cross_join(df_hazard) %>%
    #Individual hazard varies by susceptibility
    mutate(lambda = ifelse(susceptibility_status == 1,
                           lambda*hr_susceptibles,
                           lambda)) %>%
    #Indidvidual hazard might be impacted by treatment effect after treatment is initiated
    #The treatment effect might vary by susceptibility status
    mutate(hr_vaccination = ifelse(susceptibility_status == 1, hr_vaccination_susceptibles,hr_vaccination_not_susceptibles)) %>%
    mutate(lambda = ifelse(vaccine_status == 1 & time_vaccination <= time, 
                           lambda*hr_vaccination,
                           lambda)) %>%
    #Compute cumulative hazard over time - we assume hazard is piecewise constant
    group_by(id) %>%
    mutate(H = cumsum(lambda * dt))
  
  # Simulate exponential draws
  U <- runif(N)
  Z <- -log(U)
  Tevent <- numeric(N)
  
  for (i in seq_len(N)) {
    
    # find interval where cumulative hazard exceeds Z
    H_ind = c(0,data_sim$H[data_sim$id == i])
    k <- which(H_ind[-1] >= Z[i])[1]
    lambda_ind = data_sim$lambda[data_sim$id == i]
    
    if (is.na(k)) {
      # event after last knot
      Tevent[i] <- times[length(times)] +
        (Z[i] - H_ind[length(H_ind)]) / lambda_ind[length(lambda_ind)]
    } else {
      Tevent[i] <- times[k] +
        (Z[i] - H_ind[k]) / lambda_ind[k]
    }
  }
  
  # Administrative censoring
  dat_vaccination$time  <- pmin(Tevent, tau)
  dat_vaccination$event <- as.integer(Tevent <= tau)
  
  #Finally, we update vaccination status : patients who die before being vaccinated are not vaccinated anymore
  dat_vaccination = dat_vaccination %>% mutate(vaccine_status = ifelse(
    time_vaccination > time, 0,vaccine_status))
  
  return(dat_vaccination)
  
}


##############################################################
#Simulate vaccination status and vaccination times  - with piecewise constant hazard function
##############################################################

simulate_piecewise_treatment <- function(times, surv,
                                         N, tau = max(times)) {
  
  
  # Interval lengths and hazards
  dt <- diff(times)
  lambda <- log(surv[-length(surv)] / surv[-1]) / dt
  
  # cumulative hazard at knot times
  H <- c(0, cumsum(lambda * dt))
  
  # simulate exponential draws
  U <- runif(N)
  Z <- -log(U)  # total cumulative hazard
  
  Tevent <- numeric(N)
  
  for (i in seq_len(N)) {
    # find interval where cumulative hazard exceeds Z
    k <- which(H[-1] >= Z[i])[1]
    
    if (is.na(k)) {
      # event after last knot
      Tevent[i] <- times[length(times)] +
        (Z[i] - H[length(H)]) / lambda[length(lambda)]
    } else {
      Tevent[i] <- times[k] +
        (Z[i] - H[k]) / lambda[k]
    }
  }
  
  # administrative censoring
  time  <- pmin(Tevent, tau)
  event <- as.integer(Tevent <= tau)
  
  return(
    data.frame(
      id    = seq_len(N),
      time  = time,
      event = event
    ))
  
}


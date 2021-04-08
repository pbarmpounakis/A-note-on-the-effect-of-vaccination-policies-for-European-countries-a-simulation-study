#-----------------
#We will use a multitype SEIR model accounting for vaccination states for COVID-19 data in Greece
#We include 2 scenarios
#scenario 1, susceptible are given the second dose after 3 months of the first dose 
#scenario 2, susceptible are given the second dose on the recommended timeline



Stoch_multi_SEIR<-function(time_sim,parameters,initials){
  #'--------------
  #'This function performs forward simulation for the Covid-19, allowing for the implementation of Vaccination strategies
  #'We will use a stochastic multitype SEIR model, also including vaccination states (we can call it SVEIR)
  #'--------------
  #'We include 2 scenarios for the vaccinations
  #'Scenario 1: susceptibles are given the second dose after 3 periods of the first dose,
  #'having a  drop in immunity each month (length of the period and the drop in immunity are user defined)
  #'Scenario 2: susceptibles are given the second dose on the recommended time_simline
  #'
  #'
  #'-----------------
  #' 'time_sim'=length of time_sim we want to run the simulation forward
  #' -----------------
  #' 'parameters'= user defined parameters list including: transimission matrix, rates to change from different states,
  #' number of groups, population at each group,
  #' the vaccination strategies for each scenario, the number of days we assume are needed for the immunity to take effect after vaccination,
  #' the number of days we assume there is immunity drop after the first dose,  and the immunity that each group has at any given vaccination state
  #'-----------------
  #' 'initials'= vector containing the initial parameters for ech group and each state
  #'-----------------
  #'
  #'
  #'returns a dataframe with the number of people at each group and each state
  if (is.null(parameters)) {
    stop("undefined 'parameters'")
  }
  if (is.null(initials)) {
    stop("undefined 'initials'")
  }
  if (is.null(time_sim)) {
    stop("undefined 'time_sim'")
  }
  
  lambda<-parameters$transm_rates#transimission matrix
  N<-parameters$N#number of population in each group
  n_groups<-parameters$n_groups#number of groups
  vac_daily1<-parameters$vac_daily1#matrix(time_sim,n_groups) containg for each day and each group the number of people vaccinated at scenario 1
  vac_daily2<-parameters$vac_daily2#matrix(time_sim,n_groups) containg for each day and each group the number of people vaccinated at scenario 2
  days_for_immun<-parameters$days_for_immun#days for immunity after each dose
  days_immun_drop<-parameters$days_immun_drop#days for each immunity drop
  days_V22_to_V23<-parameters$days_V22_to_V23#days to 'full' immunity after the second vaccine dose in scenario 2(reqular dosing time_simline)
  immun1<-parameters$immun1/100#immunity at states V12,V22
  immun2<-parameters$immun2/100#immunity at states V13
  immun3<-parameters$immun3/100#immunity at states V14
  immun4<-parameters$immun4/100#immunity at states V15
  immun_final<-parameters$immun_final/100#immunity at states V16,V23
  #we consider fixed time_sims for E->I and I->R
  mean_time_E<-ceiling(1/parameters$rateEtoI)#mean time an individual remains at state E(integer required)
  mean_time_I<-ceiling(1/parameters$rateItoR)#mean time an individual remains at state I(integer required)
  mean_time_EV<-ceiling(1/parameters$rateEVtoIV)#mean time an individual remains at state E-Exposed but Vaccinated(integer required)
  mean_time_IV<-ceiling(1/parameters$rateIVtoRV)#mean time an individual remains at state IV-Infected but Vaccinated-less time than Infected (integer required)
  
  #define matrices for each state
  S<-matrix(as.integer(),time_sim,n_groups)#Susceptible
  E<-matrix(as.integer(),time_sim,n_groups)#Exposed
  I<-matrix(as.integer(),time_sim,n_groups)#Infected
  R<-matrix(as.integer(),time_sim,n_groups)#Removed
  EV<-matrix(as.integer(),time_sim,n_groups)#Exposed Vaccinated
  IV<-matrix(as.integer(),time_sim,n_groups)#Infected Vaccinated
  RV<-matrix(as.integer(),time_sim,n_groups)#Removed Vaccinated
  V11<-matrix(as.integer(),time_sim,n_groups)
  V12<-matrix(as.integer(),time_sim,n_groups)
  V13<-matrix(as.integer(),time_sim,n_groups)
  V14<-matrix(as.integer(),time_sim,n_groups)
  V15<-matrix(as.integer(),time_sim,n_groups)
  V16<-matrix(as.integer(),time_sim,n_groups)
  V21<-matrix(as.integer(),time_sim,n_groups)
  V22<-matrix(as.integer(),time_sim,n_groups)
  V23<-matrix(as.integer(),time_sim,n_groups)
  
  #new infections each day from every group
  new_infectionsTotal<-matrix(0,time_sim,n_groups)
  new_infectionsVac<-matrix(0,time_sim,n_groups)
  new_infectionsS<-matrix(0,time_sim,n_groups)
  new_infectionsV11<-matrix(0,time_sim,n_groups)
  new_infectionsV12<-matrix(0,time_sim,n_groups)
  new_infectionsV13<-matrix(0,time_sim,n_groups)
  new_infectionsV14<-matrix(0,time_sim,n_groups)
  new_infectionsV15<-matrix(0,time_sim,n_groups)
  new_infectionsV16<-matrix(0,time_sim,n_groups)
  new_infectionsV21<-matrix(0,time_sim,n_groups)
  new_infectionsV22<-matrix(0,time_sim,n_groups)
  new_infectionsV23<-matrix(0,time_sim,n_groups)
  #2 different  reproduction numbers
  #first eigen value of the next generation matrix
  first_eigen_value<-eigen(lambda*mean_time_I)$values[1]
  print(paste0('R0 = ',first_eigen_value))
  Reff<-rep(first_eigen_value ,time_sim)
  Rt<-rep(0,time_sim)
  # Discrete serial interval distribution
  serial_int = rep(0,time_sim)
  serial_int[1] = pgamma(1.5,shape=2.6,rate =0.4) -
    pgamma(0,shape=2.6,rate =0.4)
  for(i in 2:length(serial_int)) {
    serial_int[i] = pgamma(i+.5,shape=2.6,rate =0.4) -
      pgamma(i-.5,shape=2.6,rate =0.4)
  }
  
  
  #individuals moving from one vaccinated state to another
  V11toV12<-matrix(0,time_sim,n_groups)
  V12toV13<-matrix(0,time_sim,n_groups)
  V13toV14<-matrix(0,time_sim,n_groups)
  V14toV15<-matrix(0,time_sim,n_groups)
  V15toV16<-matrix(0,time_sim,n_groups)
  V21toV22<-matrix(0,time_sim,n_groups)
  V22toV23<-matrix(0,time_sim,n_groups)
  #initialization
  for(j in 1:n_groups){
    S[1,j]<-initials[paste0('S',j)]
    E[1,j]<-initials[paste0('E',j)]
    I[1,j]<-initials[paste0('I',j)]
    R[1,j]<-initials[paste0('R',j)]
    EV[1,j]<-initials[paste0('EV',j)]
    IV[1,j]<-initials[paste0('IV',j)]
    RV[1,j]<-initials[paste0('RV',j)]
    V11[1,j]<-initials[paste0('V11',j)]
    V12[1,j]<-initials[paste0('V12',j)]
    V13[1,j]<-initials[paste0('V13',j)]
    V14[1,j]<-initials[paste0('V14',j)]
    V15[1,j]<-initials[paste0('V15',j)]
    V16[1,j]<-initials[paste0('V16',j)]
    V21[1,j]<-initials[paste0('V21',j)]
    V22[1,j]<-initials[paste0('V22',j)]
    V23[1,j]<-initials[paste0('V23',j)]
  }
  
  for(i in 2:time_sim){
    Reff[i]<-Reff[1]*sum(S[i-1,]+V11[i-1,]+(1-immun1)*V12[i-1,]+(1-immun2)*V13[i-1,]+(1-immun3)*V14[i-1,]+(1-immun4)*V15[i-1,]+(1-immun_final)*V16[i-1,]+V21[i-1,]+(1-immun1)*V22[i-1,]+(1-immun_final)*V23[i-1,])/sum(N)
    
    for(j in 1:n_groups){
      #new infections each day
      #people in S and V11 and V21 are consider fully susceptible
      #people in other states of vaccination have some partial immunity, 
      #for this reason we multiply the probability of meeting an individual with (1-immunity)
      #R ifelse function was used to not allow for negative values in the probabilities of the binomials(it happens sometime_sims with really small values due to computer accuracy i think)
      new_infectionsS[i,j]<-rbinom(1,as.integer(S[i-1,j]),ifelse(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))>0,1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)),0))
      new_infectionsV11[i,j]<-rbinom(1,as.integer(V11[i-1,j]),ifelse(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))>0,1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)),0))
      new_infectionsV12[i,j]<-rbinom(1,as.integer(V12[i-1,j]),ifelse((1-immun1)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun1)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      new_infectionsV13[i,j]<-rbinom(1,as.integer(V13[i-1,j]),ifelse((1-immun2)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun2)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      new_infectionsV14[i,j]<-rbinom(1,as.integer(V14[i-1,j]),ifelse((1-immun3)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun3)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      new_infectionsV15[i,j]<-rbinom(1,as.integer(V15[i-1,j]),ifelse((1-immun4)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun4)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      new_infectionsV16[i,j]<-rbinom(1,as.integer(V16[i-1,j]),ifelse((1-immun_final)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun_final)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      new_infectionsV21[i,j]<-rbinom(1,as.integer(V21[i-1,j]),ifelse(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))>0,1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)),0))
      new_infectionsV22[i,j]<-rbinom(1,as.integer(V22[i-1,j]),ifelse((1-immun1)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun1)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      new_infectionsV23[i,j]<-rbinom(1,as.integer(V23[i-1,j]),ifelse((1-immun_final)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N)))>0,(1-immun_final)*(1-exp(-lambda[,j]%*%((I[i-1,]+IV[i-1,])/N))),0))
      #sum all the infections
      new_infectionsTotal[i,j]<-new_infectionsS[i,j]+new_infectionsV11[i,j]+new_infectionsV12[i,j]+new_infectionsV13[i,j]+new_infectionsV14[i,j]+new_infectionsV15[i,j]+new_infectionsV16[i,j]+new_infectionsV21[i,j]+new_infectionsV22[i,j]+new_infectionsV23[i,j]
      new_infectionsVac[i,j]<-new_infectionsV11[i,j]+new_infectionsV12[i,j]+new_infectionsV13[i,j]+new_infectionsV14[i,j]+new_infectionsV15[i,j]+new_infectionsV16[i,j]+new_infectionsV21[i,j]+new_infectionsV22[i,j]+new_infectionsV23[i,j]
      

      #we will work transmission separately from every possibly susceptible state to E
      #separately for E->I->R
      #and separately for any Vaccinated state to another Vaccinated state
      
      #-------------------------------
      #transmission portion from every possibly susceptible state to E
      
      #everyday theoretically we have vaccinations and new infections
      #so we have transmission from S->V11 and S->V21
      #and every new infection goes to E
      
      S[i,j]<-S[i-1,j]-new_infectionsS[i,j]-(vac_daily1[i,j]+vac_daily2[i,j])
      #stopping criterion if for any group the number of susceptible reaches zero
      #if the criterion is used the function will print the time_sim and the group
      #useful if we have not calculated the vaccinations properly 
      #and it avoids at the next step the function rbinom to have negative numbers,
      #because an error will be returned and the simulation won't be saved
      if(S[i,j]<=0){
        #no more vaccinations
        S[i,j]<-0
        vac_daily2[i,j]<-0
        vac_daily1[i,j]<-0
        
      }
      
      
      #---------
      #transmission portion for E->I->R not Vaccinated
      
      #1st case, nobody moves from E->I or I->R
      #2nd case we only move from E->I
      #3rd case move also from I->R
      if(i<=mean_time_E) {
        E[i,j]<-E[i-1,j]+new_infectionsS[i,j]
        I[i,j]<-I[i-1,j]
        R[i,j]<-R[i-1,j]
      }else if(i>mean_time_E & i<=(mean_time_E+mean_time_I)) {
        E[i,j]<-E[i-1,j]+new_infectionsS[i,j]- new_infectionsS[i-mean_time_E,j] 
        I[i,j]<-I[i-1,j]+new_infectionsS[i-mean_time_E,j]
        R[i,j]<-R[i-1,j]
      }else if(i>(mean_time_E+mean_time_I)){
        E[i,j]<-E[i-1,j]+new_infectionsS[i,j]- new_infectionsS[i-mean_time_E,j]
        I[i,j]<-I[i-1,j]+new_infectionsS[i-mean_time_E,j]-new_infectionsS[i-(mean_time_E+mean_time_I),j]
        R[i,j]<-R[i-1,j]+new_infectionsS[i-(mean_time_E+mean_time_I),j]
      }
      #transmission portion for EV->IV->RV  Vaccinated!!
      
      #1st case, nobody moves from EV->IV or IV->RV
      #2nd case we only move from EV->IV
      #3rd case move also from IV->RV
      if(i<=mean_time_EV) {
        EV[i,j]<-EV[i-1,j]+new_infectionsVac[i,j]
        IV[i,j]<-IV[i-1,j]
        RV[i,j]<-RV[i-1,j]
      }else if(i>mean_time_EV & i<=(mean_time_EV+mean_time_IV)) {
        EV[i,j]<-EV[i-1,j]+new_infectionsVac[i,j]- new_infectionsVac[i-mean_time_EV,j] 
        IV[i,j]<-IV[i-1,j]+new_infectionsVac[i-mean_time_EV,j]
        RV[i,j]<-RV[i-1,j]
      }else if(i>(mean_time_EV+mean_time_IV)){
        EV[i,j]<-EV[i-1,j]+new_infectionsVac[i,j]- new_infectionsVac[i-mean_time_EV,j]
        IV[i,j]<-IV[i-1,j]+new_infectionsVac[i-mean_time_EV,j]-new_infectionsVac[i-(mean_time_EV+mean_time_IV),j]
        RV[i,j]<-RV[i-1,j]+new_infectionsVac[i-(mean_time_EV+mean_time_IV),j]
      }
      #---------
      #transmission portion for the vaccinated states
      #scenario 1, people that will get the second dose after 3 periods of the first dose
      if(i<=days_for_immun){
        V11[i,j]<-V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]
        V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]
        V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
      }else if(i>days_for_immun & i<=days_immun_drop+days_for_immun){
        V11toV12[i,j]<-ifelse(vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))>0,vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j])),0)
        V11[i,j]<-ifelse(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]>0,V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j],0)
        V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]
        V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
      }else if(i>days_immun_drop+days_for_immun & i<=days_immun_drop*2+days_for_immun){
        
        
        V11toV12[i,j]<-ifelse(vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))>0,vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j])),0)
        V11[i,j]<-ifelse(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]>0,V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j],0)
        V12toV13[i,j]=ifelse(V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))>0,V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j])),0)
        V12[i,j]<-ifelse(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]>0,V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j],0)
        V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
      }else if(i>2*days_immun_drop+days_for_immun & i<=days_immun_drop*3+days_for_immun){
        
        
        V11toV12[i,j]<-ifelse(vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))>0,vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j])),0)
        V11[i,j]<-ifelse(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]>0,V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j],0)
        V12toV13[i,j]=ifelse(V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))>0,V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j])),0)
        V12[i,j]<-ifelse(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]>0,V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j],0)
        V13toV14[i,j]<-ifelse(V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j]))>0,V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j])),0)
        V13[i,j]<-ifelse(V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]>0,V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j],0)
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
      }else if(i>3*days_immun_drop+days_for_immun & i<=3*days_immun_drop+2*days_for_immun){
        
        
        V11toV12[i,j]<-ifelse(vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))>0,vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j])),0)
        V11[i,j]<-ifelse(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]>0,V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j],0)
        V12toV13[i,j]=ifelse(V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))>0,V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j])),0)
        V12[i,j]<-ifelse(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]>0,V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j],0)
        V13toV14[i,j]<-ifelse(V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j]))>0,V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j])),0)
        V13[i,j]<-ifelse(V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]>0,V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j],0)
        
        V14toV15[i,j]<-ifelse(V13toV14[i-days_immun_drop,j]-floor(mean(new_infectionsV14[(i-days_immun_drop):i,j]))>0,V13toV14[i-days_immun_drop,j]-floor(mean(new_infectionsV14[(i-days_immun_drop):i,j])),0)
        V14[i,j]<-ifelse(V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j]>0,V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j],0)
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]+V14toV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
      }else if(i>3*days_immun_drop+2*days_for_immun){
        
        
        V11toV12[i,j]<-ifelse(vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))>0,vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j])),0)
        V11[i,j]<-ifelse(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]>0,V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j],0)
        V12toV13[i,j]=ifelse(V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))>0,V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j])),0)
        V12[i,j]<-ifelse(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]>0,V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j],0)
        V13toV14[i,j]<-ifelse(V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j]))>0,V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j])),0)
        V13[i,j]<-ifelse(V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]>0,V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j],0)
        
        V14toV15[i,j]<-ifelse(V13toV14[i-days_immun_drop,j]-floor(mean(new_infectionsV14[(i-days_immun_drop):i,j]))>0,V13toV14[i-days_immun_drop,j]-floor(mean(new_infectionsV14[(i-days_immun_drop):i,j])),0)
        V14[i,j]<-ifelse(V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j]>0,V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j],0)

        V15toV16[i,j]<-ifelse(V14toV15[i-days_immun_drop,j]-floor(mean(new_infectionsV15[(i-days_immun_drop):i,j]))>0,V14toV15[i-days_immun_drop,j]-floor(mean(new_infectionsV15[(i-days_immun_drop):i,j])),0)
        V15[i,j]<-ifelse(V15[i-1,j]-new_infectionsV15[i,j]+V14toV15[i,j]-V15toV16[i,j]>0,V15[i-1,j]-new_infectionsV15[i,j]+V14toV15[i,j]-V15toV16[i,j],0)
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]+V15toV16[i,j]
      }
      #scenario 2, people that will get the second dose in a regular time_simframe
      if(i<=days_for_immun){
        V21[i,j]<-V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j]
        V22[i,j]<-V22[i-1,j]-new_infectionsV22[i,j]
        V23[i,j]<-V23[i-1,j]-new_infectionsV23[i,j]
      }
      if(i>days_for_immun & i<=days_V22_to_V23+days_for_immun){
        V21toV22[i,j]=ifelse(vac_daily2[i-days_for_immun,j]-floor(mean(new_infectionsV21[(i-days_for_immun):i,j]))>0,vac_daily2[i-days_for_immun,j]-floor(mean(new_infectionsV21[(i-days_for_immun):i,j])),0)
        V21[i,j]<-ifelse(V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i-1,j]-V21toV22[i,j]>0,V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i-1,j]-V21toV22[i,j],0)
        V22[i,j]<-V22[i-1,j]-new_infectionsV22[i,j]+V21toV22[i,j]
        V23[i,j]<-V23[i-1,j]-new_infectionsV23[i,j]
        
      } 
      if(i>days_V22_to_V23+days_for_immun){
        V21toV22[i,j]=ifelse(vac_daily2[i-days_for_immun,j]-floor(mean(new_infectionsV21[(i-days_for_immun):i,j]))>0,vac_daily2[i-days_for_immun,j]-floor(mean(new_infectionsV21[(i-days_for_immun):i,j])),0)
        V21[i,j]<-ifelse(V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i-1,j]-V21toV22[i,j]>0,V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i-1,j]-V21toV22[i,j],0)
        V22toV23[i,j]=ifelse(V21toV22[i-days_V22_to_V23,j]-floor(mean(new_infectionsV22[(i-days_V22_to_V23):i,j]))>0,V21toV22[i-days_V22_to_V23,j]-floor(mean(new_infectionsV22[(i-days_V22_to_V23):i,j])),0)
        V22[i,j]<-ifelse(V22[i-1,j]-new_infectionsV22[i,j]+V21toV22[i,j]-V22toV23[i,j]>0,V22[i-1,j]-new_infectionsV22[i,j]+V21toV22[i,j]-V22toV23[i,j],0)
        V23[i,j]<-V23[i-1,j]-new_infectionsV23[i,j]+V22toV23[i,j]
      }
      
    }
    
    #instantaneous 
    if(i==2){
      Rt[i]<-sum(new_infectionsTotal[i,])/(sum(new_infectionsTotal[1,])%*%rev(serial_int[1]))
      
    }else{
      Rt[i]<-sum(new_infectionsTotal[i,])/(rowSums(new_infectionsTotal[1:(i-1),])%*%rev(serial_int[1:(i-1)]))
    }

  }
  return(list(S=S,E=E,I=I,R=R,EV=EV,IV=IV,RV=RV,V11=V11,V12=V12,V13=V13,V14=V14,V15=V15,V16=V16,V21=V21,V22=V22,V23=V23,new_infectionsTotal=new_infectionsTotal,new_infectionsVac=new_infectionsVac,Reff=Reff,Rt=Rt))
  
}

Stoch_multi_SEIR_simulations<-function(time_sim,parameters,initials,simulations,uncertainty_interval=c(5,95)){
  #'function to perform multiple simulations(default=10) for the stochastic multitype SVEIR model
  #' same parameters as function 'Stoch_multi_SEIR' with the addition of
  #'------------
  #'  'simulations'= number of simulations to perform
  #'------------
  #'returns a dataframe just like function 'Stoch_multi_SEIR' with the means of every simulation
  
  n_groups<-parameters$n_groups#number of groups
  sim<-list()
  results<-list()
  results_median<-list()
  results_lower<-list()
  results_upper<-list()
  
  states1<-c("S","E","I","R","EV","IV","RV","V11","V12","V13","V14","V15","V16","V21","V22","V23",'new_infectionsTotal','new_infectionsVac')
  states2<-c('Reff','Rt')
  states_with_groups<-matrix(as.character(),length(states1),n_groups)
  for(j in 1:n_groups){
    states_with_groups[,j]<-paste0(states1,j)
  }
  for(state in states1){
    results_median[[state]]<-matrix(as.integer(),time_sim,n_groups)
    results_lower[[state]]<-matrix(as.integer(),time_sim,n_groups)
    results_upper[[state]]<-matrix(as.integer(),time_sim,n_groups)
    
  }
  for(state in states2){
    results_median[[state]]<-rep(0,time_sim)
    results_lower[[state]]<-rep(0,time_sim)
    results_upper[[state]]<-rep(0,time_sim)
    results[[state]]<-matrix(as.integer(),time_sim,simulations)
    
    
  }  
  for(state in states_with_groups){
    results[[state]]<-matrix(as.integer(),time_sim,simulations)
  }
  
  for(i in 1:simulations){
    sim<-Stoch_multi_SEIR(time_sim,parameters,initials)
    print(paste0('simulation ',i,' finished',sep=''))
    
    for(j in 1:length(states1)){
        for(k in 1:n_groups){
          results[[states_with_groups[j,k]]][,i]<-sim[[states1[j]]][,k]
        }
    }
    for(state in states2){
        results[[state]][,i]<-sim[[state]]
      }
    }

  
  for(j in 1:length(states1)){
    
    for(k in 1:n_groups){
      
      results_median[[states1[j]]][,k]<-apply(results[[states_with_groups[j,k]]],1,median)
      results_upper[[states1[j]]][,k]<-apply(results[[states_with_groups[j,k]]],1,quantile,probs=0.95)
      results_lower[[states1[j]]][,k]<-apply(results[[states_with_groups[j,k]]],1,quantile,probs=0.05)
      
    }
  }
  for(state in states2){
    

      results_median[[state]]<-apply(results[[state]],1,median)
      results_upper[[state]]<-apply(results[[state]],1,quantile,probs=0.95)
      results_lower[[state]]<-apply(results[[state]],1,quantile,probs=0.05)
      
    
  }
  
  return(list(results_median,results_lower,results_upper))
}

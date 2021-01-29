#-----------------
#We will use a multitype SEIR model accounting for vaccination states for COVID-19 data in Greece
#We include 2 scenarios
#scenario 1, susceptible are given the second dose after 3 months of the first dose 
#scenario 2, susceptible are given the second dose on the recommended timeline




Stoch_multi_SEIR<-function(time,parameters,initials){
  #'--------------
  #'This function performs forward simulation for the Covid-19, allowing for the implementation of Vaccination strategies
  #'We will use a stochastic multitype SEIR model, also including vaccination states (we can call it SVEIR)
  #'--------------
  #'We include 2 scenarios for the vaccinations
  #'Scenario 1: susceptibles are given the second dose after 3 periods of the first dose,
  #'having a  drop in immunity each month (length of the period and the drop in immunity are user defined)
  #'Scenario 2: susceptibles are given the second dose on the recommended timeline
  #'
  #'
  #'-----------------
  #' 'time'=length of time we want to run the simulation forward
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
  if (is.null(time)) {
    stop("undefined 'time'")
  }
  
  lambda<-parameters$transm_rates#transimission matrix
  N<-parameters$N#number of population in each group
  n_groups<-parameters$n_groups#number of groups
  vac_daily1<-parameters$vac_daily1#matrix(time,n_groups) containg for each day and each group the number of people vaccinated at scenario 1
  vac_daily2<-parameters$vac_daily2#matrix(time,n_groups) containg for each day and each group the number of people vaccinated at scenario 2
  days_for_immun<-parameters$days_for_immun#days for immunity after each dose
  days_immun_drop<-parameters$days_immun_drop#days for each immunity drop
  days_V22_to_V23<-parameters$days_V22_to_V23#days to 'full' immunity after the second vaccine dose in scenario 2(reqular dosing timeline)
  immun1<-parameters$immun1#immunity at states V12,V22
  immun2<-parameters$immun2#immunity at states V13
  immun3<-parameters$immun3#immunity at states V14
  immun4<-parameters$immun4#immunity at states V15
  immun_final<-parameters$immun_final#immunity at states V16,V23
  #we consider fixed times for E->I and I->R
  mean_time_E<-ceiling(1/parameters$rateEtoI)#mean time an individual remains at state E(integer required)
  mean_time_I<-ceiling(1/parameters$rateItoR)#mean time an individual remains at state I(integer required)
  
  #define matrices for each state
  S<-matrix(as.integer(),time,n_groups)
  E<-matrix(as.integer(),time,n_groups)
  I<-matrix(as.integer(),time,n_groups)
  R<-matrix(as.integer(),time,n_groups)
  V11<-matrix(as.integer(),time,n_groups)
  V12<-matrix(as.integer(),time,n_groups)
  V13<-matrix(as.integer(),time,n_groups)
  V14<-matrix(as.integer(),time,n_groups)
  V15<-matrix(as.integer(),time,n_groups)
  V16<-matrix(as.integer(),time,n_groups)
  V21<-matrix(as.integer(),time,n_groups)
  V22<-matrix(as.integer(),time,n_groups)
  V23<-matrix(as.integer(),time,n_groups)
  
  #new infections each day from every group
  new_infectionsTotal<-matrix(0,time,n_groups)
  new_infectionsS<-matrix(0,time,n_groups)
  new_infectionsV11<-matrix(0,time,n_groups)
  new_infectionsV12<-matrix(0,time,n_groups)
  new_infectionsV13<-matrix(0,time,n_groups)
  new_infectionsV14<-matrix(0,time,n_groups)
  new_infectionsV15<-matrix(0,time,n_groups)
  new_infectionsV16<-matrix(0,time,n_groups)
  new_infectionsV21<-matrix(0,time,n_groups)
  new_infectionsV22<-matrix(0,time,n_groups)
  new_infectionsV23<-matrix(0,time,n_groups)
  #individuals moving from one vaccinated state to another
  V11toV12<-matrix(0,time,n_groups)
  V12toV13<-matrix(0,time,n_groups)
  V13toV14<-matrix(0,time,n_groups)
  V14toV15<-matrix(0,time,n_groups)
  V15toV16<-matrix(0,time,n_groups)
  V21toV22<-matrix(0,time,n_groups)
  V22toV23<-matrix(0,time,n_groups)
  #initialization
  for(j in 1:n_groups){
    S[1,j]<-initials[paste0('S',j)]
    E[1,j]<-initials[paste0('E',j)]
    I[1,j]<-initials[paste0('I',j)]
    R[1,j]<-initials[paste0('R',j)]
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
  
  for(i in 2:time){
    for(j in 1:n_groups){
      
      #new infections each day
      #people in S and V11 and V21 are consider fully susceptible
      #people in other states of vaccination have some partial immunity, 
      #for this reason we multiply the probability of meeting an individual with (1-immunity)
      #R ifelse function was used to not allow for negative values in the probabilities of the binomials(it happens sometimes with really small values due to computer accuracy i think)
      new_infectionsS[i,j]<-rbinom(1,as.integer(S[i-1,j]),ifelse(1-exp(-lambda[j,]%*%(I[i-1,]/N))>0,1-exp(-lambda[j,]%*%(I[i-1,]/N)),0))
      new_infectionsV11[i,j]<-rbinom(1,as.integer(V11[i-1,j]),ifelse(1-exp(-lambda[j,]%*%(I[i-1,]/N))>0,1-exp(-lambda[j,]%*%(I[i-1,]/N)),0))
      new_infectionsV12[i,j]<-rbinom(1,as.integer(V12[i-1,j]),ifelse((1-immun1)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun1)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      new_infectionsV13[i,j]<-rbinom(1,as.integer(V13[i-1,j]),ifelse((1-immun2)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun2)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      new_infectionsV14[i,j]<-rbinom(1,as.integer(V14[i-1,j]),ifelse((1-immun3)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun3)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      new_infectionsV15[i,j]<-rbinom(1,as.integer(V15[i-1,j]),ifelse((1-immun4)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun4)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      new_infectionsV16[i,j]<-rbinom(1,as.integer(V16[i-1,j]),ifelse((1-immun_final)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun_final)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      new_infectionsV21[i,j]<-rbinom(1,as.integer(V21[i-1,j]),ifelse(1-exp(-lambda[j,]%*%(I[i-1,]/N))>0,1-exp(-lambda[j,]%*%(I[i-1,]/N)),0))
      new_infectionsV22[i,j]<-rbinom(1,as.integer(V22[i-1,j]),ifelse((1-immun1)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun1)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      new_infectionsV23[i,j]<-rbinom(1,as.integer(V23[i-1,j]),ifelse((1-immun_final)*(1-exp(-lambda[j,]%*%(I[i-1,]/N)))>0,(1-immun_final)*(1-exp(-lambda[j,]%*%(I[i-1,]/N))),0))
      #sum all the infections
      new_infectionsTotal[i,j]<-new_infectionsS[i,j]+new_infectionsV11[i,j]+new_infectionsV12[i,j]+new_infectionsV13[i,j]+new_infectionsV14[i,j]+new_infectionsV15[i,j]+new_infectionsV16[i,j]+new_infectionsV21[i,j]+new_infectionsV22[i,j]+new_infectionsV23[i,j]
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
      #if the criterion is used the function will print the time and the group
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
      #transmission portion for E->I->R
      
      #1st case, nobody moves from E->I or I->R
      #2nd case we only move from E->I
      #3rd case move also from I->R
      if(i<=mean_time_E) {
        E[i,j]<-E[i-1,j]+new_infectionsTotal[i,j]
        I[i,j]<-I[i-1,j]
        R[i,j]<-R[i-1,j]
      }else if(i>mean_time_E & i<=(mean_time_E+mean_time_I)) {
        E[i,j]<-E[i-1,j]+new_infectionsTotal[i,j]- new_infectionsTotal[i-mean_time_E,j] #we use E[i,j] and not E[i-1,j] because we have already add the E[i-1,j]
        I[i,j]<-I[i-1,j]+new_infectionsTotal[i-mean_time_E,j]
        R[i,j]<-R[i-1,j]
      }else if(i>(mean_time_E+mean_time_I)){
        E[i,j]<-E[i-1,j]+new_infectionsTotal[i,j]- new_infectionsTotal[i-mean_time_E,j]#we use E[i,j] and not E[i-1,j] because we have already add the E[i-1,j]
        I[i,j]<-I[i-1,j]+new_infectionsTotal[i-mean_time_E,j]-new_infectionsTotal[i-(mean_time_E+mean_time_I),j]
        R[i,j]<-R[i-1,j]+new_infectionsTotal[i-(mean_time_E+mean_time_I),j]
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
        
        
      }else if(i>days_for_immun & i<=days_immun_drop*1+days_for_immun){
        
        V11toV12[i,j]<-vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))
        if(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]<=0){
          V11[i,j]<-0
          V11toV12[i,j]<-0
        }else{
          V11[i,j]<-V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]
        }
        V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]
        V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
        
        
        
      }else if(i>days_immun_drop+days_for_immun & i<=days_immun_drop*2+days_for_immun){
        
        V11toV12[i,j]<-vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))
        if(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]<=0){
          V11[i,j]<-0
          V11toV12[i,j]<-0
        }else{
          V11[i,j]<-V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]
        }
        V12toV13[i,j]=V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))
        if(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]<=0){
          V12[i,j]<-0
          V12toV13[i,j]<-0
        }else{
          V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]
        }
        V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
        
        
        
      }else if(i>2*days_immun_drop+days_for_immun & i<=days_immun_drop*3+days_for_immun){
        
        
        V11toV12[i,j]<-vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))
        if(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]<=0){
          V11[i,j]<-0
          V11toV12[i,j]<-0
        }else{
          V11[i,j]<-V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]
        }
        V12toV13[i,j]=V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))
        if(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]<=0){
          V12[i,j]<-0
          V12toV13[i,j]<-0
        }else{
          V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]
        }
        V13toV14[i,j]<-V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j]))
        if(V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]<=0){
          V13[i,j]<-0
          V13toV14[i,j]<-0
        }else{
          V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]
        }
        V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
        
        
      }else if(i>3*days_immun_drop+days_for_immun & i<=3*days_immun_drop+2*days_for_immun){        
        
        
        
        V11toV12[i,j]<-vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))
        if(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]<=0){
          V11[i,j]<-0
          V11toV12[i,j]<-0
        }else{
          V11[i,j]<-V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]
        }
        V12toV13[i,j]=V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))
        if(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]<=0){
          V12[i,j]<-0
          V12toV13[i,j]<-0
        }else{
          V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]
        }
        V13toV14[i,j]<-V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j]))
        if(V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]<=0){
          V13[i,j]<-0
          V13toV14[i,j]<-0
        }else{
          V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]
        }
        V14toV15[i,j]<-V13toV14[i-days_immun_drop,j]-floor(mean(new_infectionsV14[(i-days_immun_drop):i,j]))
        if(V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j]<=0){
          V14[i,j]<-0
          V14toV15[i,j]<-0
        }else{
          V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j]
        }
        V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]+V14toV15[i,j]
        V16[i,j]<-V16[i-1,j]-new_infectionsV16[i,j]
        
        
        
      }else if(i>3*days_immun_drop+2*days_for_immun){
        
        
        V11toV12[i,j]<-vac_daily1[i-days_for_immun,j]-floor(mean(new_infectionsV11[(i-days_for_immun):i,j]))
        if(V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]<=0){
          V11[i,j]<-0
          V11toV12[i,j]<-0
        }else{
          V11[i,j]<-V11[i-1,j]+vac_daily1[i,j]-new_infectionsV11[i,j]-V11toV12[i,j]
        }
        V12toV13[i,j]=V11toV12[i-days_immun_drop,j]-floor(mean(new_infectionsV12[(i-days_immun_drop):i,j]))
        if(V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]<=0){
          V12[i,j]<-0
          V12toV13[i,j]<-0
        }else{
          V12[i,j]<-V12[i-1,j]-new_infectionsV12[i,j]+V11toV12[i,j]-V12toV13[i,j]
        }
        V13toV14[i,j]<-V12toV13[i-days_immun_drop,j]-floor(mean(new_infectionsV13[(i-days_immun_drop):i,j]))
        if(V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]<=0){
          V13[i,j]<-0
          V13toV14[i,j]<-0
        }else{
          V13[i,j]<-V13[i-1,j]-new_infectionsV13[i,j]+V12toV13[i,j]-V13toV14[i,j]
        }
        V14toV15[i,j]<-V13toV14[i-days_immun_drop,j]-floor(mean(new_infectionsV14[(i-days_immun_drop):i,j]))
        if(V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j]<=0){
          V14[i,j]<-0
          V14toV15[i,j]<-0
        }else{
          V14[i,j]<-V14[i-1,j]-new_infectionsV14[i,j]+V13toV14[i,j]-V14toV15[i,j]
        }
        V15toV16[i,j]<-V14toV15[i-days_immun_drop,j]-floor(mean(new_infectionsV15[(i-days_immun_drop):i,j]))
        if(V15[i-1,j]-new_infectionsV15[i,j]+V14toV15[i,j]-V15toV16[i,j]<=0){
          V15[i,j]<-0
          V15toV16[i,j]<-0
        }else{
          V15[i,j]<-V15[i-1,j]-new_infectionsV15[i,j]+V14toV15[i,j]-V15toV16[i,j]
        }
        V16[i,j]<-ifelse(V16[i-1,j]-new_infectionsV16[i,j]+V15toV16[i,j]>0,V16[i-1,j]-new_infectionsV16[i,j]+V15toV16[i,j],0)
        
        
        
      }
      #scenario 2, people that will get the second dose in a regular timeframe
      if(i<=days_for_immun){
        V21[i,j]<-ifelse(V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j]>0,V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j],0)
        V22[i,j]<-ifelse(V22[i-1,j]-new_infectionsV22[i,j]>0,V22[i-1,j]-new_infectionsV22[i,j],0)
        V23[i,j]<-ifelse(V23[i-1,j]-new_infectionsV23[i,j]>0,V23[i-1,j]-new_infectionsV23[i,j],0)
      }
      if(i>days_for_immun & i<=days_V22_to_V23+days_for_immun){
        V21toV22[i,j]=vac_daily2[i-days_for_immun,j]-floor(mean(new_infectionsV21[(i-days_for_immun):i,j]))
        if(V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j]-V21toV22[i,j]<=0){
          V21[i,j]<-0
          V21toV22[i,j]<-0
        }else{
          V21[i,j]<-V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j]-V21toV22[i,j]
        }
        V22[i,j]<-V22[i-1,j]-new_infectionsV22[i,j]+V21toV22[i,j]
        V23[i,j]<-V23[i-1,j]-new_infectionsV23[i,j]
        
      } 
      if(i>days_V22_to_V23+days_for_immun){
        V21toV22[i,j]=vac_daily2[i-days_for_immun,j]-floor(mean(new_infectionsV21[(i-days_for_immun):i,j]))
        if(V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j]-V21toV22[i,j]<=0){
          V21[i,j]<-0
          V21toV22[i,j]<-0
        }else{
          V21[i,j]<-V21[i-1,j]+vac_daily2[i,j]-new_infectionsV21[i,j]-V21toV22[i,j]
        }        
        V22toV23[i,j]=V21toV22[i-days_V22_to_V23,j]-floor(mean(new_infectionsV22[(i-days_V22_to_V23):i,j]))
        if(V22[i-1,j]-new_infectionsV22[i,j]+V21toV22[i,j]-V22toV23[i,j]<=0){
          V22[i,j]<-0
          V22toV23[i,j]<-0
          
        }else{V22[i,j]<-V22[i-1,j]-new_infectionsV22[i,j]+V21toV22[i,j]-V22toV23[i,j]
        }
        V23[i,j]<-V23[i-1,j]-new_infectionsV23[i,j]+V22toV23[i,j]
      }
      
    }
    
    
    
  }
  return(list(S=S,E=E,I=I,R=R,V11=V11,V12=V12,V13=V13,V14=V14,V15=V15,V16=V16,V21=V21,V22=V22,V23=V23,new_infectionsTotal=new_infectionsTotal))
  
}

Stoch_multi_SEIR_simulations<-function(time,parameters,initials,simulations){
  #'function to perform multiple simulations(default=10) for the stochastic multitype SVEIR model
  #' same parameters as function 'Stoch_multi_SEIR' with the addition of
  #'------------
  #'  'simulations'= number of simulations to perform
  #'------------
  #'returns a dataframe just like function 'Stoch_multi_SEIR' with the means of every simulation
  df<-list()
  df_means<-list()
  stop_time_vector<-c()
  states<-c("S","E","I","R","V11","V12","V13","V14","V15","V16","V21","V22","V23",'new_infectionsTotal')
  n_groups<-parameters$n_groups#number of groups

  for(i in 1:simulations){
    df[[i]]<-Stoch_multi_SEIR(time,parameters,initials)
    #find if any simulations finished early due to stopping criterion
  print(paste0('simulation ',i,' finished',sep=''))
  }
  
    for(j in 1:length(states)){
    sum_state<-matrix(0,time,n_groups)
for(i in 1:simulations){
  sum_state<-sum_state+df[[i]][[j]]
}

df_means[[j]]<-floor(sum_state/simulations)
  } 
    
  names(df_means)<-states
return(df_means)
}


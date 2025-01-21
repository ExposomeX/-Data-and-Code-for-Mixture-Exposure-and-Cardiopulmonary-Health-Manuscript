pacman::p_load(tidyverse,readxl,dplyr,stringr,lme4,mediation,boot,lmerTest,RMediation,data.table)

rm(list = ls())

# 1.Read ------------------------------------------------------------------

readxl::read_xlsx("data_1.OriginData.xlsx") %>% 
  as_tibble() %>% 
  dplyr::mutate(gender = as.factor(gender),
                Treat = as.factor(Treat),
                ID = as.factor(ID),
                group = as.factor(group),
                dow = as.factor(dow),
                doy = as.factor(doy),
                pollution = as.factor(pollution)) -> df_data
# 2.Tidy ------------------------------------------------------------------

df_data %>% dplyr::select(OHNAP1:DNP) %>% colnames() -> Exposure_list
df_data %>% dplyr::select(sbp:fef2575) %>% colnames() -> Outcome_list
df_data %>% dplyr::select(feno:isopgf) %>% colnames() -> Mediator_list
df_data %>% dplyr::select(gender,age,bmi,itr,t,rh) %>% colnames() -> Confusion_list

# log-transfer
for (i in Exposure_list) { 
  df_data[,i]<-log(df_data[,i],2)
}
for (i in Mediator_list) {
  df_data[,i]<-log(df_data[,i],10)
}
for (i in Outcome_list) { 
  df_data[,i]<-log(df_data[,i],10)
}

df_data %>% 
  dplyr::select(-all_of(Exposure_list)) %>% 
  cbind( 
    df_data %>% 
      dplyr::select(all_of(Exposure_list)) %>%
      mutate(across(where(is.numeric),
                    ~scales::rescale(.x)))) %>% 
  as_tibble() -> df_data_new

# 3.LME ----------------------------------------------------------

Func_LME <- function(exposure,outcome){
  formula <- as.formula(str_c(outcome,"~",str_c(c(exposure,Confusion_list),collapse = "+"),"+(1|group/ID)"))
  model <- lmerTest::lmer(formula = formula,data = df_data_new,REML = FALSE)
  a <- summary(model)
  b <- summary(model) %>% .$coefficients
  c <- confint(model,method = "Wald")
  data.frame(Exposure = exposure,Outcome = outcome,
             Coefficient = b[2,1],
             CI.Low = c[5,1], CI.High = c[5,2],
             SE = b[2,2], P.value = b[2,5]) -> temp
  
  temp %>% 
    dplyr::mutate(Coef = (exp(Coefficient)-1)*100) %>%
    dplyr::mutate(CIL = (exp(CI.Low)-1)*100) %>%
    dplyr::mutate(CIH = (exp(CI.High)-1)*100)
  return(temp)
}

expand.grid(Exposure_list,Outcome_list) %>% 
  as.data.frame() %>% 
  purrr::set_names(c("Exposure","Outcome")) -> input_vars
purrr::map2_dfr(input_vars$Exposure %>% as.vector(),
                input_vars$Outcome %>% as.vector(),Func_LME) -> df_result

df_result %>% 
  as_tibble() %>% 
  writexl::write_xlsx("Output/1.LME-Exposure&Outcome.xlsx")

# 4.Mediation ------------------------------------------------------------------
df_data_new %>% dplyr::select(all_of(Exposure_list)) -> df_exposure
df_data_new %>% dplyr::select(all_of(Outcome_list)) -> df_outcome
df_data_new %>% dplyr::select(all_of(Mediator_list)) -> df_mediator
var_exposure <- names(df_exposure)
var_outcome <- names(df_outcome)
var_mediator <- names(df_mediator)


Func_Med <- function(outcome){
  associa <- data.frame()
  for(j in 1:length(df_exposure)){
    exposure <- var_exposure[j]
    for(k in 1:length(df_mediator)){
      mediator <- var_mediator[k]
      formula_me <- as.formula(str_c(mediator,"~",str_c(c(exposure,Confusion_list),collapse = "+"),"+(1|group/ID)"))
      model_me <- lmerTest::lmer(formula_me,data = df_data_new,REML = FALSE)
      formula_out <- as.formula(str_c(outcome,"~",str_c(c(exposure,mediator,Confusion_list),collapse = "+"),"+(1|group/ID)"))
      model_out <- lmerTest::lmer(formula_out,data = df_data_new,REML = FALSE)
      
      data.frame(PredM = predict(model_me),
                 PredY = predict(model_out)) %>% 
        cbind(df_data_new) %>% 
        as_tibble() -> df
      formula_me2 <- as.formula(str_c("PredM~",str_c(c(exposure,Confusion_list),collapse = "+")))
      formula_out2 <- as.formula(str_c("PredY~",str_c(c(exposure,mediator,Confusion_list),collapse = "+")))
      
      model_me2 <- lm(formula_me2,data = df)
      model_out2 <- lm(formula_out2,data = df)
      
      mediation.result <- mediation::mediate(model.m = model_me2,model.y = model_out2,
                                             treat = var_exposure[j], mediator = var_mediator[k],
                                             boot = T, sims = 1000)
      a <- summary(mediation.result)
      associa[k+(j-1)*length(df_mediator),1] <- var_exposure[j]    
      associa[k+(j-1)*length(df_mediator),2] <- var_mediator[k]    
      associa[k+(j-1)*length(df_mediator),3] <- a$d0  
      associa[k+(j-1)*length(df_mediator),4] <- a$d0.ci[1]  
      associa[k+(j-1)*length(df_mediator),5] <- a$d0.ci[2]  
      associa[k+(j-1)*length(df_mediator),6] <- a$d0.p
      associa[k+(j-1)*length(df_mediator),7] <- a$z0
      associa[k+(j-1)*length(df_mediator),8] <- a$z0.ci[1]
      associa[k+(j-1)*length(df_mediator),9] <- a$z0.ci[2]
      associa[k+(j-1)*length(df_mediator),10] <- a$z0.p
      associa[k+(j-1)*length(df_mediator),11] <- a$tau.coef  
      associa[k+(j-1)*length(df_mediator),12] <- a$tau.ci[1]
      associa[k+(j-1)*length(df_mediator),13] <- a$tau.ci[2]
      associa[k+(j-1)*length(df_mediator),14] <- a$tau.p
      associa[k+(j-1)*length(df_mediator),15] <- a$n0  
      associa[k+(j-1)*length(df_mediator),16] <- a$n0.ci[1]
      associa[k+(j-1)*length(df_mediator),17] <- a$n0.ci[2]
      associa[k+(j-1)*length(df_mediator),18] <- a$n0.p
    }
  }
  associa %>% 
    as_tibble() %>% 
    purrr::set_names(c("Exposure","Mediator",
                       "Indirect","Indirect.CI.low","Indirect.CI.high","Indirect.p",
                       "Direct","Direct.CI.low","Direct.CI.high","Direct.p",
                       "Total","Total.CI.low","Total.CI.high","Total.p",
                       "Prop","Prop.CI.low","Prop.CI.high","Prop.p")) %>% 
    dplyr::mutate(Outcome = outcome) -> temp
  return(temp)
}
purrr::map_dfr(Outcome_list,Func_Med) %>% 
  writexl::write_xlsx("Output/2.MediationResults.xlsx")

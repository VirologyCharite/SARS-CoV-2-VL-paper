library(here)
library(boot)
library(data.table)
library(cmdstanr)
library(rjson)
library(readxl)
library(ggplot2)
library(magrittr)
library(posterior)
library(bayesplot)
library(splines)
options(mc.cores = 4)

set_cmdstan_path(
  ifelse(Sys.info()["sysname"] == "Darwin",
         "/Users/guidobiele/.cmdstanr/cmdstan-2.26.1",
         "/home/ubuntu/cmdstan-2.26.1/"))

# fulldt =
#     fread(here("data/categories-and-load-with-negatives.tsv.bz2"))
#  save(fulldt,file = here("data/fulldt.Rdata"))

#' Identify centre with longest duration between consecutive tests from the same centre
#' 
#' @param centres vector with centres
#' @param day vector with day of test (first test is day 0)
#' @return name of center with longest duration between consecutive tests from the same centre
get_ld_centre = function(centres,day) {
  streaks = data.frame(centre = rep("",10), duration = rep(NA,10))
  k = 0
  duration = 0
  for (d in 2:length(centres)) {
    if (centres[d] == centres[d-1] & d < length(centres)) {
      duration = duration + day[d] - day[d-1]
    } else if (centres[d] != centres[d-1]) {
      k = k+1
      streaks$centre[k] = centres[d-1]
      streaks$duration[k] = duration
      duration = 0
    } else if (centres[d] == centres[d-1] & d == length(centres)) {
      duration = duration + day[d] - day[d-1]
      k = k+1
      streaks$centre[k] = centres[d-1]
      streaks$duration[k] = duration
    }
  }
  streaks = streaks[complete.cases(streaks),]
  streaks = streaks[order(-streaks$duration),]
  
  if (max(streaks$duration) == 0) {
    tbl = sort(table(streaks$centre),decreasing = T)
    return(
      ifelse(tbl[1] != tbl[2],
             names(tbl)[1],"X")
    )
  } else if (nrow(streaks) == 1 | streaks$duration[1] != streaks$duration[2]) {
    return(streaks$centre[streaks$duration == max(streaks$duration)])
  } else {
    tmp = paste0(streaks[streaks$duration == max(streaks$duration),"centre"], collapse = "_")
    if (grepl("ICU",tmp)) {
      return("ICU")
    } else if (grepl("IDW",tmp)) {
      return("IDW") 
    } else if ( (!grepl("IDW",tmp)) & (!grepl("ICU",tmp)) & grepl("WD",tmp)) {
      return("WD")
    } else if (grepl("C19",tmp) & grepl("ED",tmp)) {
      return("ED")
    }
  }
}

grp_var = "ID"
HOSPITAL_TEST_CENTRE = c('H', 'ICU', 'IDW', 'WD')
PAMS_TEST_CENTRE = c('C19')

#### Load & prepare data for analysis ####

#' Load an preprocess full data set
#' 
#' @return a data.table 
get_full_data = function() {
  AG_breaks = c(seq(0,25,5),seq(35,65,10),101)
  AG_labels = c(paste0(c(seq(0,20,5),seq(25,55,10)),"-",c(seq(5,25,5),seq(35,65,10))),">65")
  my_breaks = c(seq(0,25,5),seq(35,65,10),120)
  my_labels = c("0-5","5-10","10-15","15-20","20-25","25-35","35-45","45-55","55-65",">65")
  load(here("data/fulldt.Rdata"))
  fulldt = 
    fulldt %>% 
    .[Age <= 0, Age := 0.1] %>%
    .[, AgeGroup := cut(Age,breaks = AG_breaks,ordered = T, labels = AG_labels)] %>%
    .[, Date := as.Date(Date,format = "%Y-%m-%d")] %>%
    .[, month := month(Date)] %>%
    .[, day := 1+as.numeric(Date - min(Date))] %>%
    .[, PCR := gsub("T2","cobas",PCR)] %>%
    .[, PCR := factor(PCR)] %>%
    setnames("PCR1-TestCentreCategory","TestCentreCategory") %>%
    .[, PCR := factor(PCR)] %>%
    .[ TestCentreCategory == "C19" & Hospitalized == 0, PAMS1 := 1] %>% 
    .[ , Hospital_centre := ifelse(TestCentreCategory %in% HOSPITAL_TEST_CENTRE,1,0)] %>% 
    .[, Group := factor(ifelse(PAMS1 == 1, "PAMS",ifelse(Hospitalized == 1,"Hospitalized","Other")),
                        levels = c("Other","PAMS","Hospitalized"))] %>%
    .[, PAMS := ifelse(TestCentreCategory %in% PAMS_TEST_CENTRE, 1, 0)] %>% 
    .[, `Age category` := cut(Age, breaks = my_breaks, ordered_result = T,
                              labels = my_labels)] %>% 
    setkeyv(c("Age","personHash"))
  
  AgeGrpLvls = fulldt[, list(paste(floor(range(Age)), collapse = "-")), by = "AgeGroup"]$V1
  fulldt %>% 
    .[, Month := ordered(months(Date),
                         levels = unique(months(sort(fulldt$Date))))]
  
  fulldt %>% 
    .[, nAge := Age + rnorm(1, sd = .5), by = 1:nrow(fulldt)] %>% 
    .[, B117centre := ifelse(TestCentreFull %in% unique(fulldt[B117 == 1, TestCentreFull]),T,F)] %>% 
    .[, cWeek := week(Date)] %>% 
    .[, CentreWeek := paste0(TestCentreFull,cWeek)]
  
  return(fulldt)
}

#' Load an preprocess first positive test data
#' 
#' @return a data.table 
get_log10Load_data = function() {
  bdata = 
    get_full_data() %>% 
    .[log10Load > 0]
  
  
  B117CentreDay = 
    unique(bdata[B117 == 1, .(TestCentreFull,Date)])
  
  bdata %>% 
    .[, B117CentreDay0 := F] %>% 
    .[, B117CentreDay1 := F] %>% 
    .[, B117CentreDay2 := F] %>% 
    .[, B117CentreDay3 := F] %>% 
    .[, B117CentreDay4 := F] %>% 
    .[, B117CentreDay5 := F]
  for (k in 1:nrow(B117CentreDay)) {
    for (d in 0:5) {
      bdata[TestCentreFull == B117CentreDay$TestCentreFull[k] & 
              abs(as.numeric(Date-B117CentreDay$Date[k])) <= d,
            (paste0("B117CentreDay",d)) := T]
    }
  }
  
  contrasts(bdata$AgeGroup,2) = contr.poly(length(unique(bdata$AgeGroup)))
  
  return(bdata)
}

#' Load an preprocess culture positivity data
#' 
#' @return a data.table 
get_CP_data = function() {
  dt = 
    read_excel(here("data/2020-08-19_merged_extracted_data.xlsx")) %>%
    data.table() %>%
    setnames(c("day_from_symptom_onset_raw",   # renaming variables
               "cycle_threshold_extracted",
               "log10_viral_load_extracted",
               "original_study"),
             c("day_from_onset",
               "Ct",
               "log10Load",
               "Study"))
  
  wdt = rbind(
    fread(here("data/2020-08-27_woelfel_isolation_data.csv")) %>%
      .[, log10_viral_load := as.numeric(gsub(",",".",log10_viral_load, perl = T))] %>%
      .[, Study := "woelfel"] %>%
      .[, Ct := NA] %>%
      setnames(c("log10_viral_load","days_since_symptom","patient"),
               c("log10Load","day_from_onset","ID"))
    ,
    fread(here("data/2020-08-20_woelfel_swab_sputum_dirty_wide.csv")) %>%
      melt(id.var = "days_since_symtpom", value.name = "log10Load") %>%
      .[, variable := as.character(variable)] %>%
      .[!is.na(log10Load)] %>%
      .[, ID := strsplit(variable, "_")[[1]][1],
        by = c("days_since_symtpom","variable","log10Load")] %>%
      .[, sample_type := ifelse(grepl("swap",variable),"swab","sputum")] %>%
      .[, variable := NULL] %>%
      .[, Study := "woelfel"] %>% 
      .[, culture_outcome := NA] %>%
      setnames(c("days_since_symtpom"),
               c("day_from_onset"))
    , fill = T) %>%
    .[, Ct := log10Load * -3.962 + 49.678]
  
  setkeyv(wdt,c("ID","sample_type","day_from_onset"))
  
  
  dt = 
    rbind(dt[,intersect(names(dt),names(wdt)), with = F],
          wdt[,intersect(names(dt),names(wdt)), with = F]) %>%
    .[culture_outcome == "no_culture", culture_positive := NA] %>%
    .[, culture_positive := ifelse(culture_outcome == "positive",1L,0L)] 
  
  ber_dt = fromJSON(file = here("data/sa-min-3.json"))
  ber_dt = do.call(rbind,lapply(ber_dt$people,
                                function(x)
                                  as.data.table(x))) %>%
    .[, ID := factor(personHash)] %>%
    .[, ID := factor(as.numeric(ID))] %>%
    .[, c("firstName","familyName") := NULL] %>%
    .[, PCRdate := as.Date(PCRdate, format = "%m/%d/%Y")]  %>%
    .[, day := as.numeric(PCRdate)] %>% 
    .[, DOB := as.Date(DOB, format = "%m/%d/%Y")]  %>%
    .[, Age := as.numeric((PCRdate-DOB))/365] %>%
    .[, Age := min(Age), by = "ID"] %>%
    setnames("viralLoad", "log10Load") %>%
    .[, N_tests := .N, by = ID] %>%
    .[, Study := "BER"]

  for (id in unique(ber_dt$ID)) 
    ber_dt[ID == id, day := day - min(ber_dt[ID == id,day])]

  load(here("data/fulldt.Rdata"))
  setkeyv(fulldt,"personHash")
  setkeyv(ber_dt,"personHash")
  ber_dt[fulldt, `:=`(B117 = B117, Gender = Gender, PAMS1 = PAMS1, Hospitalized = Hospitalized)]
  
  merged_data =
    rbind(dt,ber_dt, fill= T) %>%
    .[!is.na(log10Load)] %>%
    .[, Study := factor(Study)] %>%
    .[, StudyID := as.numeric(Study)] %>% 
    .[, N_tests := .N, by = "ID"] %>%
    .[, N_test_grp := cut(N_tests,breaks = c(0,2.5,3,4,5,6,7,9,15,100), ordered_result = T)] %>%
    .[, N_test_grp := cut(N_tests,breaks = c(0,5,15,100), ordered_result = T)] %>%
    .[,c("ID","day_from_onset", "Study", "log10Load", "culture_positive","N_test_grp","N_tests","Age","day","testName",
         "centreCategory","Gender","Hospitalized","PCRdate","personHash","B117")] 
  return(merged_data)
}

#' Preprocess data for viral load time course analysis
#' 
#' @return a data.table 
prep_time_course_data = function(merged_data, latest_peak_day = Inf) {
  PAMS1_personHash = get_log10Load_data()[PAMS1 == 1,personHash]
  ## cleaning options
  # 1. remove time courses with peak load after day 15
  
  tmp = 
    merged_data %>%
    .[,PAMS1 := 0] %>%
    .[personHash %in% PAMS1_personHash, PAMS1 := 1] %>%
    .[, Gender := as.numeric(factor(Gender,levels = c("F","M")))-1] %>% 
    .[Study == "BER"] %>%
    .[, N_tests := .N, by = ID] %>%
    .[ N_tests >= 3] %>%
    setkeyv(unique(c(grp_var,"ID","day_from_onset"))) %>%
    .[, max_day := max(day), by = "ID"] %>%
    .[, max_load := max(log10Load), by = "ID"] %>%
    .[, max_load_day := max((log10Load == max_load)*day), by = .(ID)] %>% 
    .[log10Load > 2, day_last_positive := max(day), by = ID] %>%
    .[, day_last_positive := max(day_last_positive,na.rm = T), by = ID] %>%
    .[log10Load > 2, day_first_positive := min(day), by = ID] %>%
    .[, day_first_positive := max(day_first_positive,na.rm = T), by = ID] %>%
    .[, days := max_day - min(day_first_positive, na.rm = T), by = "ID"] %>%
    .[, hospitalized := sum(centreCategory %in% HOSPITAL_TEST_CENTRE) > 0, by = "ID"] %>%
    .[, phosptests := mean(centreCategory %in% HOSPITAL_TEST_CENTRE), by = "ID"] %>%
    .[, NtestsCat := cut(N_tests, breaks = c(2:9,20))] %>%
    .[, last_test_negative := ifelse(tail(log10Load,1) <= 2 ,T,F), by = "ID"] %>%
    .[, first_test_negative := ifelse(head(log10Load,1) <= 2 ,T,F), by = "ID"] %>%
    .[, first_last_test_negative := ifelse(last_test_negative == T & first_test_negative == T,T,F), by = "ID"] %>% 
    .[, ld_centre := get_ld_centre(centreCategory,day), by = "ID"] %>% 
    .[, diff_load_12 := diff(log10Load)[1], by = .(ID)] %>% 
    .[, diff_day_12 := diff(day)[1], by = .(ID)] %>% 
    .[, diff_load_12perday := diff_load_12/diff_day_12, by = .(ID)] %>% 
    .[, Test1_positive_increasing := ifelse(diff_load_12perday > 0 & first_test_negative == F,1,0)] %>% 
    .[, Test1_positive_decreasing := ifelse(diff_load_12perday <= 0 & first_test_negative == F,1,0)] %>% 
    .[, Test1_positive_malo := ifelse(diff_load_12perday < 0 & first_test_negative == F,max_load,0)] %>% 
    .[, Test1_pos_incr_malo := ifelse(diff_load_12perday > 0 & first_test_negative == F,max_load,0)] %>% 
    .[, Test1_pos_decr_malo := ifelse(diff_load_12perday < 0 & first_test_negative == F,max_load,0)] %>% 
    .[, Test1_negative_malo := ifelse(first_test_negative == T,max_load,0)] %>% 
    .[max_load_day < latest_peak_day]
    
  tmp[ld_centre %in% names(which(table(tmp[day == 0, ld_centre])<21)), ld_centre := "X"]
  
  
  return(tmp)
}

#' Generate list with data for Stan model for time course analysis
#' 
#' @param selection minimum number of tests per subject
#' @param samples reduced number of subjects chose (for testing)
#' @param max_diff_load_12perday maximum initial raw slope allowed. Will deselect participants with a higher slope
#' @param ub_log_slope_up_mu upper bound for the Stan model parameter log_slope_up_mu
#' @return a list
make_time_course_standata = function(selection = 3,
                                     samples = NULL,
                                     max_diff_load_12perday = NULL,
                                     ub_log_slope_up_mu = Inf,
                                     imputation_limit = 3,
                                     latest_peak_day = Inf) {
  merged_data = get_CP_data()
  day_data = 
    prep_time_course_data(merged_data, latest_peak_day = latest_peak_day)
  
  
  if (!is.null(max_diff_load_12perday)) {
    exclude_ID = day_data[diff_load_12perday > max_diff_load_12perday & first_test_negative == T & day == 0,ID] 
    if (length(exclude_ID) > 0) {
      exclude_idx = which(day_data$ID %in% exclude_ID & day_data$day == 0)
      day_data = day_data[-exclude_idx]
      day_data %>% 
        .[, N_tests := .N, by = .(ID)] %>% 
        .[, min_day := min(day), by = .(ID)] %>% 
        .[, day := day-min_day]
    }
  }
  
  day_data = 
    day_data %>%
    #.[ ! (ID %in% unique(day_data[day > 25 & log10Load > 7,ID]))] %>% 
    .[ N_tests >= as.numeric(selection)]  
  
  if (!is.null(samples))
    day_data = day_data[ID %in% sample(day_data[day == 0,ID],samples)]
  
  lvls = names(which(table(day_data[, get(grp_var)]) > 0))
  ld_centre_lvls = c("WD",sort(setdiff(unique(day_data$ld_centre),"WD")))
  day_data %>% 
    .[, (grp_var) :=  factor(get(grp_var), levels = lvls)] %>% 
    .[, ld_centre := factor(ld_centre, levels = ld_centre_lvls)]
  
  
  centre_matrix = model.matrix(~ 0 + centreCategory, day_data[day == 0])
  Gender = day_data[day == 0, Gender]
  X_PGH = model.matrix(~ 1 + PAMS1 + Gender + Hospitalized + B117 + PAMS1:Hospitalized + B117:PAMS1, day_data[day == 0][, Gender := ifelse(is.na(Gender),0.5,Gender)])
  X_PGH = X_PGH[,-1]
  X_PGH = X_PGH[,colSums(X_PGH)>0]
  B_Age = t(bs(day_data[day == 0, Age],
               degree=3,
               knots = quantile(day_data[day == 0, Age],
                                probs = seq(.05,.95, length.out = 5),
                                names = F)))
  B_phosptests = t(bs(day_data[day == 0, phosptests],
                      degree=3,
                      knots = quantile(day_data[day == 0, phosptests],
                                       probs = seq(.05,.95, length.out = 3),
                                       names = F)))
  X_ld_centre = model.matrix(~ ld_centre, day_data[day == 0])[,-1]

  FTN = 
  datalist_DAY = list(
    N_DAY = nrow(day_data),
    gstart_DAY = c(1,which(diff(as.numeric(day_data[, get(grp_var)])) != 0)+1),
    gend_DAY = c(which(diff(as.numeric(day_data[, get(grp_var)])) != 0),nrow(day_data)),
    G = length(unique(day_data[,get(grp_var)])),
    Y_DAY = day_data$log10Load,
    X_DAY = day_data$day,
    N_T1_neg = sum(day_data[day == 0,first_test_negative] == T),
    T1_neg_idx = which(day_data[day == 0,first_test_negative] == T),
    N_T1_pos = sum(day_data[day == 0,first_test_negative] == F),
    T1_pos_idx = which(day_data[day == 0,first_test_negative] == F),
    DAY_max_load = day_data[day == 0, max_load_day],
    N_NegTests = sum(day_data$log10Load == 0),
    idx_NegTests = which(day_data$log10Load == 0),
    PCR = 1*(day_data$testName == "T2"),
    N_centres = length(unique(day_data$centreCategory)),
    centre = model.matrix(~ 0 + centreCategory, day_data),
    Gender = do.call(c,lapply(Gender, function(x) ifelse(is.na(x),0.5,1*(x == "M")))),
    N_centre1 = length(unique(day_data[day == 0, centreCategory])),
    centre1 = as.numeric(factor(day_data[day == 0, centreCategory])),
    N_ld_centre = length(unique(day_data[day == 0, ld_centre])),
    ld_centre = as.numeric(factor(day_data[day == 0, ld_centre])),
    K_ld_centre = ncol(X_ld_centre),
    X_ld_centre = X_ld_centre,
    imputation_limit = imputation_limit,
    num_basis_Age = nrow(B_Age),
    B_Age = B_Age,
    Age = as.numeric(scale(day_data[day == 0, Age])),
    num_basis_phosptests = nrow(B_phosptests),
    B_phosptests = B_phosptests,
    phosptests = day_data[day == 0, phosptests],
    X_PG = X_PGH[,1:2],
    K_PGH = ncol(X_PGH),
    X_PGH = X_PGH,
    max_load = as.numeric(scale(day_data[day == 0, max_load])),
    ub_log_slope_up_mu = ub_log_slope_up_mu
  )
  
  #### CP data ###
  cp_data = 
    merged_data[!is.na(culture_positive) & Study %in% c("woelfel","ranawaka") & log10Load > 2] 
  
  datalist_CP = list(
    N_CP = nrow(cp_data),
    Y_CP = cp_data$culture_positive,
    X_CP = matrix(cp_data$log10Load,ncol = 1)
  )
  
  ### Day and CP analysis ###
  datalist = c(datalist_CP,datalist_DAY)
  datalist$condition_on_data = 1
  
  return(list(datalist = datalist, day_data = day_data))
} 

make_age_group_N = function(bdata,breaks = NULL,my_group = NULL) {
  complete_age_table = 
    bdata[,.(Age)] %>% 
    .[, Age := floor(Age)] %>%
    .[, age_group := cut(Age, breaks = breaks, right = F)] %>% 
    unique()
  
  tmpdata = bdata[,.(Age,Group)]
  if (!is.null(my_group))
    tmpdata = tmpdata[Group %in% my_group]
  
  age_group_N =
    tmpdata %>%
    .[, Age := floor(Age)] %>%
    .[, age_group := cut(Age, breaks = breaks, right = F)] %>%
    .[, .(N_Age = sum(.N)), by = .(age_group,Age)] 
  
  setkeyv(complete_age_table,names(complete_age_table))
  setkeyv(age_group_N,names(complete_age_table))
  age_group_N = 
    complete_age_table[age_group_N, N_Age := N_Age] %>% 
    .[is.na(N_Age), N_Age := 0]
  return(age_group_N)
}

#### Manipulate draws and calculate statistics ####

#' Generate a draws data.table
#' 
#' @param draws a draws object from the posterior package
#' @param thin optional parameter to select subset of posterior draws
#' @return a data.table
as_draws_dt = function(draws, thin = 1) {
  dt = 
    as_draws_df(draws) %>%
    thin_draws(thin) %>%
    data.table() %>% 
    .[,c(".chain",".iteration") := NULL]
  return(dt)
}

#' Generate list with basic stats (mean 5% and 95% quantiles). Mainly for use in data.tables
#' 
#' @param x a vector with real numbers.
#' @return a list
post_stats_list = function(x) {
  m = collapse::fmean(x)
  qs = quantile(x,c(.05,.95), names = F)
  return(list(m = m, q5 = qs[1], q95 = qs[2]))
}

#' Generate list with mean and a larger number of quantiles. Mainly for use in data.tables
#' 
#' @param x a vector with real numbers.
#' @param quantiles a vector with quantiles to be calculates (default : seq(.5,.95,by = .05))
#' @return a list
my_stats_list_long = function(x, quantiles = seq(.5,.90,by = .05)) {
  m = collapse::fmean(x)
  qs = quantile(x,c((1-quantiles)/2,1-(1-quantiles)/2), names = F)
  dstats = as.list(c(m = m, qs))
  names(dstats) = c("mean",c(paste0("lower",quantiles*100),paste0("upper",quantiles*100)))
  return(dstats)
}

#' Calculate posterior statistics for a number of variables
#' 
#' @param dt a long-form data.table with posterior samples and auxiliary variables
#' @param var variable in data.frame for which statistics are calculated
#' @param by grouping variable in data.frame
#' @param quantiles a vector with quantiles to be calculates (default : seq(.5,.95,by = .05))
#' @return a data.table
get_stats = function(dt, var = "value", by = NULL,  quantiles = seq(.5,.95,by = .05))  {
  my_stats = 
    dt[, as.list(my_stats_list_long(get(var),quantiles = quantiles)),
       by = by]
  return(my_stats)
}

#' Extract posterior draws and group by person identifiers (id)
#' 
#' @param draws a draws object from the posterior package
#' @param params parameter to be extracted from the draws object
#' @param thin 
#' @return a data.table
draws_by_id = function(draws,params, thin = 1) {
  for (p in params) {
    tmp = 
      subset_draws(draws,p) %>%
      as_draws_dt(thin = thin) %>%
      melt(id.vars = ".draw", variable.name = "ID") %>%
      setnames("value",p)
    
    if (length(params) == 1 | p == params[1]) {
      draws.dt = copy(tmp)
    } else {
      draws.dt = cbind(draws.dt, tmp[,c(p),with = F])
    }
  }
  draws.dt %>% 
    .[, ID := as.numeric(gsub("[^0-9]","",ID,perl = T))] %>%
    .[, ID := factor(ID, labels = levels(day_data[, ID]))] %>%
    setkeyv(c(".draw","ID"))
  return(draws.dt)
}

#' Calculate statics over posterior draws that were first aggregated
#' 
#' @param dt a data.table with an indicator variable `.draw`
#' @param by grouping variable
#' @param target.var parameter for which statistics are calculated
#' @return a data.table
summarise_draws_dt_by = function(dt,by, target.var = "value", varname = NULL) {
  summarised_draws = 
    dt %>%
    .[, list(m = collapse::fmean(get(target.var))), by = c(".draw",by)] %>%
    .[, as.list(post_stats_list(m)),
      by = by]
  if (!is.null(varname))
    setnames(summarised_draws,"m",varname)
  return(summarised_draws)
}

#' Summarise draws by person identifier by using the summarise_draws function from the posterior package.
#' 
#' @param draws a draws object from the posterior package
#' @param var parameter to be extracted
#' @return a data.table with ID, parameter name and statistics from the posterior
smrs_by_ID = function(draws, var) {
  tmp = 
    subset_draws(draws,var) %>% 
    summarise_draws() %>%
    data.table() %>%
    .[, ID := as.numeric(gsub("[^0-9]","",variable, perl = T))] %>%
    .[, ID := factor(ID, labels = levels(day_data[, ID]))] %>%
    .[, parameter := var] %>%
    .[, c("ID","parameter","mean","q5","q95","sd")]
}

#' Summarise draws by a grouping variable
#' 
#' @param draws a draws object from the posterior package
#' @param var parameter to be extracted
#' @param grp.dt a data table with a column "ID" and an additional column with a grouping variable
#' @return a data.table with statistics from the posterior
smrs_by_grp = function(draws, var, grp.dt, calc.delta = F) {
  grp.var = setdiff(names(grp.dt),"ID")
  tmp = 
    draws_by_id(draws, var)
  setkeyv(tmp,"ID")
  setkeyv(grp.dt,"ID")
  tmp = 
    tmp[grp.dt, c(grp.var) := get(grp.var)] %>%
    .[, list(value = collapse::fmean(get(var))), by = c(".draw",grp.var)] 
  tmp[[grp.var]] = as.character(tmp[[grp.var]])
  if (calc.delta == T) {
    delta = tmp %>% 
      dcast(as.formula(paste(".draw ~ ", grp.var)),
            value.var = "value")
    g1 = names(delta)[2]
    g2 = names(delta)[3]
    delta %>% 
      .[, value := get(g1) - get(g2)] %>% 
      .[, (grp.var) := paste0(g1,"-",g2)]
    tmp = rbind(tmp,
                delta[,c(".draw",grp.var,"value"),with = F])
  }
  tmp %>%  
    get_stats(by = grp.var) %>%
    .[, parameter := var] %>% 
    .[,tbl := paste0(round(mean,ifelse(parameter == "slope_down",3,2)),
                     " (", round(lower90,ifelse(parameter == "slope_down",3,2)),", ",
                     round(upper90,ifelse(parameter == "slope_down",3,2)),")")] 
}

#' Calculate mean and sd over individual lecel parameters for each draw
#' 
#' @param draws a draws object from the posterior package
#' @param var parameter to be extracted. Has to be a parameter with id-specific values.
#' @return a data.table with statistics (mean, sd over IDs) from the posterior for each draw.
smrs_by_draw = function(draws, var) {
  tmp = 
    subset_draws(draws,var) %>% 
    as_draws_dt() %>%
    melt(id.vars = ".draw") %>%
    .[, c("param","ID") := tstrsplit(variable,"\\[")] %>%
    .[, ID := as.numeric(gsub("[^0-9]","",ID, perl = T))] %>%
    .[, variable := NULL] %>%
    .[, list(mean = collapse::fmean(value),
             sd = collapse::fsd(value)), by = ".draw"] %>%
    .[, parameter := var]
}

#' Calculate statistics over posteriors
#' 
#' @param draws a draws object from the posterior package
#' @param file a file with posterior draws or cmdstan model object.
#' @param vars variable for which statistics are calculated
#' @return a data.table with statistics (mean, quantiles) from the posterior.
stats_over_draws = function(draws = NULL, file = NULL, vars) {
  if (is.null(draws)) {
    load(file)
    if (is.null(draws)) draws = csf$draws()
  }
  return(
    do.call(rbind,
            lapply(vars, function(x) {
              subset_draws(draws,x) %>% 
                as_draws_dt() %>%
                melt(id.vars = ".draw") %>%
                .[, variable := NULL] %>%
                .[, list(value = collapse::fmean(value)), by = ".draw"] %>%
                get_stats() %>%
                .[, parameter := x]
            }))
  )
}


#' Generate predicted time courses of viral load and culture positivity for each individual and posterior sample
#' 
#' @param draws a draws object from the time course model generated with the posterior package
#' @param days vector of days (can be fractions of days) for which to calculate predictions
#' @param thin 
#' @return a data.table with time courses of viral load and culture positivity
make_VLCP_by_draw_ID = function(draws, days = seq(-10,30,by = 1), thin = 1) {
  by_draw_ID = 
    draws_by_id(draws,
                c("slope_up","slope_down","intercept"),
                thin = thin)
  
  beta_sweight_draw = 
    subset_draws(draws,c("beta_sweight_mu")) %>% 
    as_draws_dt(thin = thin) 
  
  by_draw_ID = merge(by_draw_ID,
                     beta_sweight_draw,
                     by = c(".draw"),
                     allow.cartesian = T)
  
  CP_params = subset_draws(draws,c("alpha_CP","beta_CP")) %>%
    as_draws_dt(thin = thin)  %>%
    setnames("beta_CP[1]","beta_CP")
  
  yhat = function(intercept, slope_up, slope_down, beta) {
    wght_d = inv.logit(days*beta)
    return((intercept + days*slope_up) * (1-wght_d) + (intercept + days*slope_down) * wght_d)
  }
  
  
  VLCP_by_draw_ID = 
    by_draw_ID[,
               as.list(yhat(intercept,slope_up,slope_down,beta_sweight_mu)),
               by = c("ID",".draw")] %>%
    melt(id.var = c("ID",".draw"), value.name = "log10Load") %>%
    .[, variable := as.numeric(gsub("[^0-9]","",variable))] %>%
    .[, day_shifted := days[variable]] %>%
    .[, variable := NULL]
  
  setkeyv(VLCP_by_draw_ID,".draw")
  setkeyv(CP_params,".draw")
  VLCP_by_draw_ID = 
    VLCP_by_draw_ID[CP_params, CP := inv.logit(alpha_CP + log10Load * beta_CP)] 
  
  return(VLCP_by_draw_ID)
}

#' Calculate posterior predictions for first positive test
#' 
#' @param brms fit object
#' @param ndt new data with columns PCR_Group TestCentreCategory Gender PCR Group and Age-Group
#' @param wghts weights to weight according to Age
#' @return a data.table with posterior predictions by Age and PAMS
calc_post_lin_pred = function(bfit,ndt,wghts, sum.Group = "Group") {
  setkeyv(wghts,"TestCentre.Group")
  ndt = ndt[TestCentre.Group %in% unique(wghts$TestCentre.Group)]
  group.by =  c("Group","Age",".draw")
  if (is.null(sum.Group)) group.by =  c("Age",".draw")
  post_lin_pred = 
    ndt %>% 
    cbind(t(posterior_epred(bfit,
                            newdata = ndt))) %>%
    melt(id.vars = names(ndt),
         variable.name = ".draw",
         value.name = "log10Load") %>%
    .[, .draw := as.numeric(gsub("V","",.draw, perl = T))] %>%
    setkeyv("TestCentre.Group") %>% 
    .[wghts, weighted_load := log10Load * weight ] %>%
    .[, .(log10Load = sum(weighted_load)),
      by = group.by]
}

#' Calculate percent positive culture by age while correctly weighting individual for each age group
#' 
#' @param post_lin_pred data.table with posterior linear predictions of log10 viral load
#' @param weights data.table with weights 
#' @param w.var name of variable with weights 
#' @return a data.table with the posterior distribution of culture positivities by age
predict_VLCP_by_Age = 
  function(post_lin_pred,weights,w.var, get.stats = T) {
    stats_by_age = 
      post_lin_pred %>% 
      .[weights, `:=`(wVL = log10Load * get(w.var), wCP = CP * get(w.var))] %>% 
      .[, .(log10Load = sum(wVL), CP = sum(wCP)), by = .(Age,.draw)] %>% 
      melt(id.vars = c("Age",".draw"), variable.name = "outcome") 
    
    if (get.stats == T) {
      return(
        stats_by_age %>% 
          get_stats(var = "value", by = c("Age", "outcome")))
    } else {
      return(stats_by_age)
    }
  }

#' Calculate viral load or culture positivity difference
#' 
#' @param dt a data.table with posterior predictions for ages 0:100
#' @param age_group_N data.table with number of subjects per age year
#' @param CPstat indicator to calculate either risk differences (RD) or risk rations (RR)
#' @return a data.table statistics for all differences
calc_diff = function(dt, age_group_N, CPstat = "RD") {
  age_group_N %>% 
    .[, N_age_group := sum(N_Age), by = .(age_group)] %>% 
    .[, w := N_Age / N_age_group]
  setkeyv(age_group_N,"Age")
  setkeyv(dt,"Age")
  dt = 
    dt %>%
    .[age_group_N, w_value := value * w] %>%
    .[age_group_N, age_group := age_group] %>%
    .[, .(value = sum(w_value)),
      by = .(age_group,outcome,.draw)] %>%
    dcast(outcome + .draw ~ age_group)
  
  if (length(levels(age_group_N$age_group)) == 5) {
    dt[, `0-5 vs 20-100` := `[0,5)` - `[20,101)`]  %>%
      .[, `5-10 vs 20-100` := `[5,10)` - `[20,101)`]  %>%
      .[, `10-15 vs 20-100` := `[10,15)` -  `[20,101)`] %>%
      .[, `15-20 vs 20-100` := `[15,20)` - `[20,101)`]
  } else if (length(levels(age_group_N$age_group)) == 6) {
    dt[, `0-5 vs 20-65` := `[0,5)` - `[20,65)`]  %>%
      .[, `5-10 vs 20-65` := `[5,10)` - `[20,65)`]  %>%
      .[, `10-15 vs 20-65` := `[10,15)` -  `[20,65)`] %>%
      .[, `15-20 vs 20-65` := `[15,20)` - `[20,65)`]
  }
  else {
    dt[, `0-5 vs 45-55` := `[0,5)` - `[45,55)`]  %>%
      .[, `5-10 vs 45-55` := `[5,10)` - `[45,55)`]  %>%
      .[, `10-15 vs 45-55` := `[10,15)` -  `[45,55)`] %>%
      .[, `15-20 vs 45-55` := `[15,20)` - `[45,55)`] %>%
      .[, `20-25 vs 45-55` := `[20,25)` - `[45,55)`] %>%
      .[, `25-35 vs 45-55` := `[25,35)` - `[45,55)`] %>%
      .[, `55-65 vs 45-55` := `[55,65)` - `[45,55)`] %>%
      .[, `65-80 vs 45-55` := `[65,80)` - `[45,55)`] %>%
      .[, `80-100 vs 45-55` := `[80,101)` - `[45,55)`]
  }
  
  dt %>% melt(id.vars = c("outcome",".draw"),
              variable.name = "comparison")
}

#' Get parameter names to calculate conditional effects for time  course model
#' 
#' @param var variable for which conditional effects are calculated
#' @param pred predictor for which conditional effects are calculated
#' @return a list with parameter names
get_cond_effect_params = function(var,pred) {
  intercept_var = paste0("log_",var,"_mu")
  a0_var = paste0("a0_",var,"_",pred)
  a_var = paste0("a_",var,"_",pred)
  b_var = paste0("betaPGH_",var)
  
  intercept = 
    draws %>% 
    subset_draws(intercept_var) %>% 
    as_draws_matrix()
  a0 = 
    draws %>% 
    subset_draws(a0_var) %>% 
    as_draws_matrix()
  a = draws %>% 
    subset_draws(a_var) %>% 
    as_draws_matrix()
  b = 
    draws %>% 
    subset_draws(b_var) %>% 
    as_draws_matrix()
  
  return(list(intercept = intercept,
              a0 = a0,
              a = a,
              b = b))
}

#' Calculate conditional effects for time course model for a continuous predictors
#' 
#' @param var variable for which conditional effects are calculated
#' @param pred predictor for which conditional effects are calculated
#' @return `data.table` with mean and quantiles for conditional predictions
cond_eff_cont = function(var, predictor = "Age") {
  
  ppost = get_cond_effect_params(var,predictor)
  
  X_PGH = datalist[["X_PG"]]
  X_PG[,1] = mean(X_PG[,1])
  X_PG[,2] = mean(X_PG[,2])
  B_predictor = datalist[[paste0("B","_",predictor)]]
  
  yhat = 
    cbind(
      day_data[day == 0, c(predictor,"ID"), with = F],
      do.call(rbind,
              lapply(1:dim(ppost$intercept)[1], 
                     function(j) {
                       ppost$intercept[j] + 
                         ppost$a0[j] * datalist[[predictor]] + 
                         ppost$a[j,] %*%  B_predictor + 
                         t(X_PG %*% t(ppost$b[j,])) }
              )
      ) %>% 
        t()
    ) %>%
    melt(id.vars = c(predictor,"ID"), variable.name = ".draw")
  
  if (predictor == "Age")
    yhat[, Age := round(Age,1)]
  yhat = 
    yhat[, list(value = collapse::fmean(value)), by = c(".draw",predictor)]
  
  yhat[, value := exp(value)]
  if (var == "slope_down")
    yhat[, value := -value]
  
  yhat[, parameter := var]
  return(yhat)
}

#' Calculate conditional effects for time course model for a continuous predictors
#' 
#' @param var variable for which conditional effects are calculated
#' @param pred predictor for which conditional effects are calculated
#' @return `data.table` with mean and quantiles for conditional predictions
cond_eff_contPGH = function(var, predictor = "Age") {
  
  ppost = get_cond_effect_params(var,predictor)
  
  X_PGH = datalist[["X_PGH"]]
  for (k in 1:ncol(X_PGH)) 
    X_PGH[,k] = mean(X_PGH[,k])
  B_predictor = datalist[[paste0("B","_",predictor)]]
  my_predictor = datalist[[predictor]]
  
  yhat = 
    cbind(
      day_data[day == 0, c(predictor,"ID"), with = F],
      do.call(rbind,
              lapply(1:dim(ppost$intercept)[1], 
                     function(j) {
                       ppost$intercept[j] + 
                         ppost$a0[j] * my_predictor + 
                         ppost$a[j,] %*%  B_predictor + 
                         t(X_PGH %*% t(ppost$b[j,])) }
              )
      ) %>% 
        t()
    ) %>%
    melt(id.vars = c(predictor,"ID"), variable.name = ".draw")
  
  if (predictor == "Age")
    yhat[, Age := round(Age,1)]
  yhat = 
    yhat[, list(value = collapse::fmean(value)), by = c(".draw",predictor)]
  
  yhat[, value := exp(value)]
  if (var == "slope_down")
    yhat[, value := -value]
  
  yhat[, parameter := var]
  return(yhat)
}

#' Calculate conditional effects for time course model for binary predictors
#' 
#' @param var variable for which conditional effects are calculated
#' @param pred predictor for which conditional effects are calculated
#' @return `data.table` with mean and quantiles for conditional predictions
cond_eff_bin = function(var, predictor = "PAMS1", Age_conditioned = mean(day_data[day == 0, Age])) {
  ppost = get_cond_effect_params(var,"Age")
  
  X_PG = datalist[["X_PG"]] 
  X_PG[,setdiff(colnames(X_PG),predictor)] = mean(X_PG[,setdiff(colnames(X_PG),predictor)])
  X_PG = unique(X_PG)
  B_Age = datalist[["B_Age"]]
  Age = datalist[["Age"]]
  idx = which.min(abs(day_data[day == 0, Age] - Age_conditioned))[1]
  B_Age = matrix(B_Age[,idx],ncol = 1)
  Age = Age[idx]
  
  yhat = 
    cbind(
      data.table(X_PG)[,predictor, with = F],
      do.call(rbind,
              lapply(1:dim(ppost$intercept)[1], 
                     function(j) {
                       as.numeric(ppost$intercept[j,]) + 
                         as.numeric(ppost$a0[j,] * Age) + 
                         as.numeric(ppost$a[j,] %*% B_Age) + 
                         t(X_PG %*% t(ppost$b[j,])) }
              )
      ) %>% 
        t()
    ) %>%
    melt(id.vars = predictor, variable.name = ".draw")
  
  yhat[, value := exp(value)]
  if (var == "slope_down")
    yhat[, value := -value]
  
  yhat[, parameter := var]
  return(yhat)
}

#### ggplot formatting ####

theme_set(
  theme_minimal() + 
    theme(axis.ticks.length = unit(2,"pt"),
          axis.title = element_text(size = 11,
                                    margin = margin(0,0,0,0)),
          axis.text = element_text(size = 10,
                                   margin = margin(0,0,0,0)),
          strip.text.x = element_text(size = 11,
                                      margin = margin(0,0,.5,0)),
          plot.background = element_rect(linetype = "blank"),
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ))

my_clrs = c("black","#DC0000FF","#4DBBD5FF")
my_clrs = c("black","red","blue")

#' Add manual scales with 3 colors, black, rad, and blue, for fill and color to ggplot plot
#' 
#' @return `list` list with manual scales for fill and color
black_red_blue = function() {
  list(
    scale_colour_manual(values = my_clrs), 
    scale_fill_manual(values = my_clrs))
}

#' Add manual scales with two colors red and blue, for fill and color to ggplot plot
#' 
#' @return `list` list with manual scales for fill and color
red_blue = function(o = 2:3) {
  list(
    scale_colour_manual(values = my_clrs[o]), 
    scale_fill_manual(values = my_clrs[o]))
}

#' Add manual scales with two colors, green and brown,  for fill and color to ggplot plot
#' 
#' @return `list` list with manual scales for fill and color
green_brown = function() {
  list(
    scale_colour_manual(values = c("#00A087FF","#7E6148FF")), 
    scale_fill_manual(values = c("#00A087FF","#7E6148FF")))
}

#' Adjust text size in ggplots
#' 
#' @return `ggplot-theme` with modified font sizes for various text elements
gg_text_size = function(scale_size = 1) {
  theme(axis.text=
          element_text(
            size=8 * scale_size,#5.5,
            margin = c(0,0,0,0, units = "pt")),
        axis.title = 
          element_text(
            size=8 * scale_size,#6,
            margin = c(0,-5,0,0, units = "pt")),
        strip.text.x = element_text(size=6 * scale_size),#6), 
        strip.text.y = element_text(size=6 * scale_size),#6), 
        plot.title = element_text(size = 10 * scale_size,#7,
                                  face = "bold")) 
  return()
}

#' Adjust legend size in ggplots
#' 
#' @return `list` with ggplot-theme for modified legend text and icons
gg_legend_size = function(ncol = 3) {
  return(
    list(
      theme(legend.text = element_text(size = 8),#5),
            legend.title = element_text(size = 9),#5),
            legend.key.size = unit(0.5, "cm")),
      guides(fill = guide_legend(ncol = ncol,
                                 keywidth = .75,
                                 keyheight = .75),
             color = guide_legend(ncol = ncol,
                                  keywidth = .75,
                                  keyheight = .75))))
}

#' Remove axes from ggplot
#' 
#' @return `theme` with axis elements set to element_blank()
theme_marginal = function() {
  theme(axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.background = 
          element_rect(fill = 
                         adjustcolor("black",alpha = 0)))
}

#' Expand x and y axes in ggplot
#' 
#' @return `list` with instructions to expand axes
gg_expand = function() {
  list(
    scale_x_continuous(expand = expansion(0,0)), 
    scale_y_continuous(expand = expansion(.01,0)))
}

#' Add grid lines to ggplot
#' 
#' @param axis is a character X, y, or xy indicating for which axis grid lines are added
#' @return `theme` with instructions to add grod lines
gg_add_grid = function(axis = "xy") {
  if (axis == "xy") {
    my_grid = theme(panel.grid.major = element_line(colour="grey", size=0.05))
  } else if (axis == "x") {
    my_grid = theme(panel.grid.major.x = element_line(colour="grey", size=0.05))
  } else if (axis == "y") {
    my_grid = theme(panel.grid.major.y = element_line(colour="grey", size=0.05))
  } 
  return(my_grid)
}

#' Add progressively shaded credible interval lines to ggplot
#' 
#' @param data is a `data.table` or `data.frame` in wide format with mean and quantiles 
#' lower50, lower55, lower95 ,... upper50, upper55, upper95 in columns 
#' @param size is the line width
#' @param alpha is the transparency level for each ribbon
#' @return a ggplot layer with a progressively shaded confidence ribbon 
conf_linerange = function(data, size = 1.5, alpha = .05) {
  lapply(seq(50,90,5), 
         function(x) 
           geom_linerange(alpha = alpha,
                          aes_string(ymin = paste0("lower",x), ymax = paste0("upper",x)),
                          position = position_dodge(width = dodge_with),
                          size = size))
} 

#' Add progressively shaded credible interal ribbons to ggplot
#' 
#' @param data is a `data.table` or `data.frame` in wide format with mean and quantiles 
#' lower50, lower55, lower95 ,... upper50, upper55, upper95 in columns 
#' @param fill is a character indicating the fill color or a variable in the data
#' @param alpha is the transparency level for each ribbon
#' @return a ggplot layer with a progressively shaded confidence ribbon 
conf_ribbon = function(data, fill = "red", alpha = .05) {
  conf_levels = grep("lower",names(data), value = T) %>% gsub("lower","",.)
  if (exists(fill,data)) {
    lapply(conf_levels, 
           function(x) 
             geom_ribbon(alpha = alpha, color = NA,
                         aes_string(ymin = paste0("lower",x), ymax = paste0("upper",x),fill = fill)))
  } else {
    lapply(conf_levels, 
           function(x) 
             geom_ribbon(alpha = alpha, fill = fill,color = NA,
                         aes_string(ymin = paste0("lower",x), ymax = paste0("upper",x))))
  }
} 

#### Plotting results ####


#' Plot time course of viral load or culture positivity.
#' each individual is plotted as a blue line and 
#' group level mean and credible intervals are plotted in red
#' 
#' @param data is a `data.table` or `data.frame` in wide format with mean and quantiles 
#' lower50, lower55, lower95 ,... upper50, upper55, upper95 in columns 
#' @param by_day is a data.table with group level statistics
#' @param by_day_id is a data.table with individual level means
#' @param y.var is the name of the variable on the y-axis
#' @param clr is the color for group level results
#' @return a ggplot layer with a progressively shaded confidence ribbon 
plot_by_day = function(by_day, by_day_id = NULL, y.var, xlim = c(-7.5,27.5), ylim = NULL , clr = "red") {
  if (exists(clr,by_day)) {
    p = ggplot(by_day, aes_string(x = "day_shifted", y = y.var, color = clr)) 
    gl = geom_line()
  } else {
    p = ggplot(by_day, aes_string(x ="day_shifted", y = y.var)) 
    gl = geom_line(color = "red")
  }
  
  if (!is.null(by_day_id)) {
    p = 
      p + geom_line(data = by_day_id, aes_string(group = "ID"), alpha = .075, color = "blue", size = .15)
  }
  
  p = 
    p + 
    coord_cartesian(xlim = xlim, ylim = ylim) + 
    gl + 
    conf_ribbon(by_day, fill = clr) +
    xlab("Days from peak viral load") + 
    gg_expand()
  
  return(p)
}

#' Plot histograms of posterior distributions 
#' 
#' @param drws is an object (`posterior` `draws` or `data.table`)
#' @param labels data.table with text to be shown on histogram (thought statistics)
#' @param fill variable name that determines fill color. if `is.null( fill)` the fill color is blue 
#' @param value.var is the name of the variable for which histogram is drawn
#' @param nrow is the number of rows in the `facet_wrap`
#' @return a ggplot plot
plot_post_hists = function(drws, labels = NULL, fill = NULL, value.var = "value", nrow = 1,labeller = label_value) {
  
  if (!("data.table" %in% class(drws))) {
    drws = 
      drws %>% 
      as_draws_dt() %>%
      melt(id.var = ".draw") 
  }
  
  if (!is.null(labels))
    drws[, variable := factor(variable, labels = labels)]
  
  if (is.null(fill)) {
    p = 
      drws %>%
      ggplot(aes_string(x = value.var)) + 
      geom_histogram(bins = 30, fill = adjustcolor("blue", alpha = .5)) 
  } else {
    p = 
      drws %>%
      ggplot(aes_string(x = value.var, fill = fill)) + 
      geom_histogram(bins = 30, alpha = .5) +
      red_blue()
  }
  
  if (length(unique(drws$variable)) > 0)
    p = p + facet_wrap(~variable, scales = "free", nrow = nrow, labeller = labeller) 
  
  gb = ggplot_build(p)
  
  label_data = 
    gb$data[[1]] %>% 
    data.table() %>%
    .[, list(x = min(x), y = max(y)), by = "PANEL"] %>%
    .[, variable := factor(PANEL, labels = sort(unique(drws$variable)))]
  
  tmp_stats = 
    drws[, as.list(post_stats_list(get(value.var))),
         by = variable] %>%
    .[, label := paste0(round(m,ifelse(m>1,1,2)),
                        " (",round(q5,ifelse(m>1,1,2)),", ",
                        round(q95,ifelse(m>1,1,2)),")")] 
  
  label_data = 
    merge(label_data[,c("variable","x","y")],
          tmp_stats[,c("variable","label")],
          by = "variable")
  
  p = 
    p + 
    theme(axis.line.y = element_blank(), axis.text.y = element_blank()) + 
    geom_text(data = label_data,hjust = 0,
              aes(x = x, y = y, label = label, fill = NULL),
              size = 4) + 
    ylab("Posterior draws")
  
  return(p)
}

#' Plot 2-d posterior predictive plot for first positive analysis
#' 
#' @param fulldt are the observed data
#' @param yhat is the posterior predictions from model
#' @param value.var is the name of the variable for which histogram is drawn
#' @param strat_var is a grouping variable
#' @return a faceted plot of 2d posterior prediction distributions
ppc_2d = function(fulldt, yhat, strat_var){
  dt = 
    cbind(data.table(fulldt[,get(strat_var)]) %>% 
            setnames("V1",strat_var),
          t(yhat)) %>%
    melt(id.vars = strat_var,
         variable.name = "iter",
         value.name = "log10Load") %>%
    .[, list(mean = mean(log10Load),
             sd = collapse::fsd(log10Load)),
      by = c(strat_var,"iter")]
  
  obs = 
    fulldt %>% 
    .[, list(mean = mean(log10Load),
             sd = collapse::fsd(log10Load)),
      by = strat_var]
  
  # specify desired contour levels:
  prob <- c(0.95,0.90,.8,0.5)
  n = 100
  
  nrows = n^2*length(unique(fulldt[,get(strat_var)]))
  dc = data.table(mean = rep(0.0,nrows),
                  sd = rep(0.0,nrows),
                  value = rep(0.0,nrows),
                  prob = rep(0.0,nrows)) %>%
    .[, (strat_var) := fulldt[1,get(strat_var)]]
  k = 0
  for (g in unique(fulldt[,get(strat_var)])) {
    mv.kde <- kde2d(dt[get(strat_var) == g, mean], dt[get(strat_var) == g,sd], n = n)
    dx <- diff(mv.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
    dy <- diff(mv.kde$y[1:2])
    sz <- sort(mv.kde$z)
    c1 <- cumsum(sz) * dx * dy
    
    dimnames(mv.kde$z) <- list(mv.kde$x,mv.kde$y)
    tmp <- reshape::melt(mv.kde$z)
    names(tmp)[1:2] = c("mean","sd")
    idx = k*n^2+(1:n^2)
    dc[idx, mean := tmp$mean]
    dc[idx, sd := tmp$sd]
    dc[idx, value := tmp$value]
    dc[idx, prob := approx(sz,1-c1,tmp$value)$y]
    dc[idx, (strat_var) := g]
    k = k+1
  }
  
  facet_formula = as.formula(paste0("~",strat_var))
  
  p = 
    ggplot(dc,aes(x=sd,y=mean))+
    geom_contour(aes(z=prob,color=..level..),breaks=prob) +
    geom_hline(data = obs, aes(yintercept = mean),
               col = adjustcolor("red",alpha = .25)) +
    geom_vline(data = obs, aes(xintercept = sd),
               col = adjustcolor("red",alpha = .25)) +
    geom_point(data=obs,col = "red",size=1) +
    facet_wrap(facet_formula) + 
    labs(color = "CI level")
  
  return(p)
}

#' Plot group differences in viral load or culture positivity of time course analysis
#' 
#' @param by_day_draws is a `data.table` with time courses for each individual and posterior draw
#' @param y.var the variable for which the difference is calculated 
#' @param grp.dt is a data.table with a columns with IDs and a second column with the grouping variable
#' @param stat is "RD" and "RR" for risk difference and risk ration, respectively 
#' @return a ggplot with group differences and credible intervals
plot_delta_grps = function(by_day_draws, y.var = "log10Load", grp.dt, stat = "RD", comparisons = NULL, plot = T) {
  grp.var = setdiff(names(grp.dt),"ID")
  cast_formula = as.formula(paste0(".draw ~", grp.var))
  
  setkeyv(by_day_draws,"ID")
  setkeyv(grp.dt,"ID")
  tmp = 
    by_day_draws[day_shifted == 0] %>% 
    .[grp.dt, c(grp.var) := get(grp.var)] %>%
    .[,.(value = mean(get(y.var))), by = c(".draw",grp.var)] %>%
    .[!is.na(get(grp.var))] %>%
    dcast(cast_formula)
  if (is.null(comparisons)) {
    comparisons = combn(names(tmp)[-1],2)
  }
  
  comp_stats = c()
  for (k in 1:ncol(comparisons)) {
    if (stat == "RR") { ## risk ratio
      delta = tmp[, get(comparisons[1,k])] / tmp[, get(comparisons[2,k])]  
    } else if (stat == "RD") { ## risk difference
      delta = tmp[, get(comparisons[1,k])] - tmp[, get(comparisons[2,k])]  
    }
    comp_stats = rbind(comp_stats,
                       get_stats(data.table(delta = delta), var = "delta"))
  }
  rm(tmp)
  tmp = gc()
  comp_stats = 
    data.table(comp_stats) %>%
    .[, comparison := paste0(comparisons[1,],
                             ifelse(stat == "RD","–","/"),
                             ifelse(ncol(comparisons) > 22,"\n",""),comparisons[2,])] %>%
    .[, x := ncol(comparisons):1]
  
  if (stat == "RR") {
    ylim = c(0.5,1.1)
  } else  if (stat == "RD") {
    offset = (max(comp_stats[,upper90])-min(comp_stats[,lower90]))*.5
    ylim = c(min(comp_stats[,lower90])-offset,max(comp_stats[,upper90]))
  }
  
  if (plot == T) {
    return(
      ggplot(comp_stats, aes(x = x, y = mean)) + 
        geom_hline(yintercept = ifelse(stat == "RR",1,0),
                   col = "red", lty = 3, size = .5) + 
        geom_point(color = "black") + 
        conf_linerange() + 
        theme(axis.line.y = element_blank(),
              axis.text.y = element_blank()) + 
        geom_text(aes(y = min(comp_stats[,lower90]) - offset*.75,
                      label = comparison),
                  size = 2.5,
                  hjust = 0)  +
        xlab("Comparison") + 
        ylab(stat) + 
        coord_cartesian(ylim = ylim) + 
        coord_flip()
    )
  } else {
    return(
      comp_stats
    )
  }
  
}

#' Plot time courses of viral load or culture positivity by group
#' 
#' @param by_day_draws is a `data.table` with time courses for each individual and posterior draw
#' @param y.var the variable for which the difference is calculated 
#' @param grp.dt is a data.table with a columns with IDs and a second column with the grouping variable
#' @param stat is "RD" and "RR" for risk difference and risk ration, respectively 
#' @return a ggplot with time courses and credible intervals
plot_by_day_grp = function(by_day_grp, grp.var, grp.dt, y.var = "mean") {
  
  p = plot_by_day(by_day_grp,
                  clr = grp.var,
                  y.var = y.var) 
  if (length(unique(grp.dt[[grp.var]])) == 2) {
    p = p + red_blue()
  }
  
  N_table = 
    table(grp.dt[,get(grp.var)]) %>%
    data.table() %>%
    .[ N > 0] %>%
    .[, label := paste0(V1," (",N,")")]
  
  p = 
    p + 
    theme(legend.position = c(.8,.8)) + 
    gg_legend_size(1)
  
  legend_title = paste(gsub("_"," ",grp.var), " (N)")
  if (class(grp.dt[,get(grp.var)])[1] == "ordered") {
    p = p +
      scale_color_ordinal(name = legend_title, labels = N_table$label) + 
      scale_fill_ordinal(name = legend_title, labels = N_table$label)
  } else {
    if (nrow(N_table) == 2) {
      color_values = c("red","blue")
    } else if (nrow(N_table) == 3) {
      color_values = c("black","red","blue")
    } else {
      color_values = colorspace::qualitative_hcl(nrow(N_table), palette = "Dark 3")
    }
    p = p +
      scale_color_manual(name = legend_title, labels = N_table$label, values = color_values) + 
      scale_fill_manual(name = legend_title, labels = N_table$label, values = color_values) 
  }
  return(p)
}

#' Combined plot of time courses of viral load and culture positivity by group and group differences 
#' 
#' @param by_day_draws is a `data.table` with time courses for each individual and posterior draw
#' @param grp.dt is a data.table with a columns with IDs and a second column with the grouping variable
#' @param stat is "RD" and "RR" for risk difference and risk ration, respectively 
#' @param comparisons optional matrix with pairwise comparisons
#' @return a ggplot with three panels: viral load time course, culture positivity time course and group differences
plot_by_day_grp_delta = function(by_day_draws,grp.dt, stat = "RD", comparisons = NULL) {
  grp.var = grep("ID",names(grp.dt), value = T, invert = T)
  by_day_grp = 
    do.call(rbind,
            lapply(unique(grp.dt[[grp.var]]), 
                   function(x) {
                     idx = grp.dt[get(grp.var) == x,ID]
                     by_day_draws[ID %in% idx,
                                  list(log10Load = collapse::fmean(log10Load),
                                       CP = collapse::fmean(CP)),
                                  by = .(.draw,day_shifted)] %>% 
                       melt(id.vars = c(".draw", "day_shifted")) %>% 
                       get_stats(by = c("variable","day_shifted")) %>% 
                       .[,(grp.var) := x]
                   })
    )
  
  VL = plot_by_day_grp(
    by_day_grp[variable == "log10Load"], grp.var = grp.var, grp.dt = grp.dt) + 
    ylab(expression(log[10]~viral~load)) +
    coord_cartesian(ylim = c(2,9), xlim = c(-8,25)) +
    theme(legend.position = c(.75,.825))
  tmp = gc()
  CP = plot_by_day_grp(
    by_day_grp[variable == "CP"], grp.var = grp.var, grp.dt = grp.dt) + 
    ylab("Probability of positive culture") + 
    coord_cartesian(xlim = c(-8,25), ylim = c(-.005,1)) +
    theme(legend.position = c(.75,.825))
  tmp = gc()
  Peak_delta = 
    plot_delta_grps(by_day_draws[day_shifted == 0],
                    grp.dt = grp.dt,
                    y.var = "CP",
                    stat = stat,
                    comparisons = comparisons) + 
    ylab("Difference of peak culture probability")
  tmp = gc()
  return(VL + CP + Peak_delta)
}

#' Plot prior and posterior of a parameter from the time course model
#' 
#' @param parameter is the name of the parameter
#' @param draws are the posterior draws from the time course model
#' @param mu is the mean of the prior distribution 
#' @param sigma is the standard deviation of the prior distribution 
#' @param prior_dist is prior distribution, which is either "normal" (default) or "log-normal"
#' @return a ggplot with the prior distribution as a black pdf and the posterior distribution as a blue density
plot_prior_posterior = function(parameter, draws, mu, sigma, prior_dist = "normal", xlim = NULL) {
  
  title = paste0("Prior: ",
                 ifelse(prior_dist == "normal","N",prior_dist),
                 "(",mu,",",sigma,")")
  
  p = 
    draws %>% 
    subset_draws(parameter) %>% 
    as_draws_dt() %>%
    setnames(parameter,"x") %>% 
    ggplot(aes(x = x)) + 
    geom_density(fill = adjustcolor("blue", alpha = .5)) +
    ggtitle(title) + 
    gg_expand() + 
    xlab(parameter) +
    theme(plot.title  = element_text(size = 10))
  
  if (is.null(xlim))
    xlim = subset_draws(draws,parameter)
  my_data = data.frame(x = seq(xlim[1],xlim[2],length.out = 101))
  if (prior_dist == "normal") {
    p = p +
      stat_function(data = my_data,
                    fun = dnorm,
                    args = list(mean = mu, sd = sigma))
  } else if (prior_dist == "log-normal") {
    p = p +
      stat_function(data = my_data,
                    fun = dlnorm,
                    args = list(meanlog = mu, sdlog = sigma))
  }
  return(p)
}

#' Plots viral load or culture positivity for a specific number of days
#' 
#' @param by_dayIDdraw data.table with draws by ID and day
#' @param target.var variable to be plotted
#' @param grp.var grouping variable
#' @param key_days days for which viral load or culture positivity are plotted
#' @return a ggplot with the viral load or culture positivity for selected days
plot_key_days_by_age = function(by_dayIDdraw, target.var = "value", var, grp.var = "Age_group", key_days = c(-2, 0, 5, 10), sub_sample_ids = NULL) {
  if (is.null(sub_sample_ids))
    sub_sample_ids = unique(by_dayIDdraw$ID)
  
  grp.dt = unique(day_data[day == 0,c("ID",grp.var),with = F])
  
  if( !exists(grp.var,by_dayIDdraw)) {
    by_dayIDdraw = merge(
      by_dayIDdraw,
      grp.dt, by = "ID")
  }
  
  by_dayAgeGrp = 
    by_dayIDdraw[day_shifted %in% key_days & ID %in% sub_sample_ids] %>%
    .[, list(value = collapse::fmean(get(target.var))), by = c("day_shifted",grp.var,".draw")] %>%
    get_stats(by = c("day_shifted",grp.var)) %>%
    .[, days_after_peak_load := ordered(day_shifted)] %>%
    setnames("mean",var)
  
  plot_by_Age = 
    ggplot(by_dayAgeGrp, aes_string(x = "days_after_peak_load", y = var, color = grp.var, group = grp.var)) + 
    geom_point(position = position_dodge(width = dodge_with), size = .5) + 
    conf_linerange(by_dayAgeGrp) + 
    xlab("Days from peak viral load") +
    theme(legend.position = c(.7,1),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4)) +
    guides(color = guide_legend(ncol = 2,
                                title = "Age group",
                                keywidth = .25,
                                keyheight = .25))
  
  N_table = 
    table(grp.dt[,get(grp.var)]) %>%
    data.table() %>%
    .[ N > 0] %>%
    .[, label := paste0(V1," (",N,")")]
  
  plot_by_Age = 
    plot_by_Age + 
    theme(legend.position = c(.8,.8)) + 
    gg_legend_size(2)
  
  legend_title = paste(gsub("_"," ",grp.var), " (N)")
  if (class(grp.dt[,get(grp.var)])[1] == "ordered") {
    plot_by_Age = plot_by_Age +
      scale_color_ordinal(name = legend_title, labels = N_table$label) + 
      scale_fill_ordinal(name = legend_title, labels = N_table$label)
  } else {
    plot_by_Age = plot_by_Age +
      scale_color_discrete(name = legend_title, labels = N_table$label) + 
      scale_fill_discrete(name = legend_title, labels = N_table$label)
  }
  
  gc()
  return(plot_by_Age)
}


#### Misc ####
sprint_stat = function(draws, digits = 1, my_quantiles = c(.05,.95)) {
  m = mean(draws)
  qs = quantile(draws,my_quantiles)
  return(sprintf(paste0("%.",digits,"f (%.",digits,"f, %.",digits,"f)"),
                 m,
                 qs[1],
                 qs[2]))
}

sprint_mci = function(x, digits = 1) {
  m = mean(x)
  cis = c(m-se(x)*1.96,m+se(x)*1.96)
  return(sprintf(paste0("%.",digits,"f (%.",digits,"f, %.",digits,"f)"),
                 m,
                 cis[1],
                 cis[2]))
}

sprint_mcib = function(x, digits = 1) {
  m = mean(x)
  se = sqrt((m/(1-m))/length(x))
  cis = c(m-se*1.96,m+se*1.96)
  return(sprintf(paste0("%.",digits,"f (%.",digits,"f, %.",digits,"f)"),
                 m*100,
                 max(cis[1]*100,0),
                 min(cis[2]*100,100)))
}

print_gstats= function(s,o,cx) {
  age_group_diff_tbl[sample == s & outcome == o][[cx]] 
}
clean_colnames = function(x) {
  colnames(x) = gsub("~ ","",colnames(x)) %>% gsub("\\[|\\]","~",.)
  return(x)
}

se = function(x) {return(sd(x)/sqrt(length(x)-1))}

#' Generate initial values to fit Stan model for age-viral load association
#' 
#' @param bfit brms fit object of model that shall be fitted again. this is here only used to get parameter names and dimensions.
#' @return a list with inital values
make_age_fit_inits = function(bfit) {
  pars = grep("^r_|^lp__|^s_|_Intercept$",names(bfit$fit@par_dims),invert = T, value = T)
  inits = vector(length = length(pars),
                 mode = "list")
  names(inits) = pars
  for (p in pars) {
    pdim = bfit$fit@par_dims[[p]]
    if (length(pdim) == 0) {
      if (grepl("^imp_|^sd_|^sds_|^Intercept$",p)) {
        inits[[p]] = runif(1)/10
      } else {
        inits[[p]] = rnorm(1,0,.1)
      }
    } else {
      if (grepl("imp_gender",p)) {
        inits[[p]] = array(runif(cumprod(pdim)),dim = pdim)
      } else if (grepl("^sd_|^sds_",p)) {
        inits[[p]] = array(runif(cumprod(pdim)),dim = pdim)/10
      } else {
        inits[[p]] = array(rnorm(cumprod(pdim),0,.1),dim = pdim)
      }
    }
  }
  inits[["Intercept"]] = runif(1,5,7)
  if (exists("nu",inits)) inits[["nu"]] = runif(length(inits[["nu"]]), 3, 20)
  return(inits)
}

#' Get statistics from brms model for B.1.1.7 analysis
#' 
#' @param bfit brms fit object 
#' @param window time window in which non-B.1.1.7 cases need to be to be compared to B.1.1.7 cases
#' @param model string with name of the model
#' @param paired indicator if only centres are used that reported both B.1.1.7 and non-B.1.1.7 cases
#' @return data table with results and description of model and sample
get_B117_stats = function(bfit,model,window,paired = "No") {
  effect.B117 = 
    bfit$fit %>% 
    as_draws() %>% 
    subset_draws("b_B117", regex = T) %>% 
    as.numeric() %>% 
    sprint_stat(2)
  
  return(
    data.table(
      Window = window,
      Adjusted = ifelse(grepl("adjusted|re.adjusted",model),"Yes","No"),
      `Rand. eff`  = ifelse(grepl("re",model),"Yes","No"),
      Paired =  paired,
      `N B.1.1.7` = sum(bfit$data$B117 == "B117" | bfit$data$B117 == 1),
      `N non-B.1.1.7` = nrow(bfit$data) - sum(bfit$data$B117 == "B117" | bfit$data$B117 == 1),
      `Effect B.1.1.7` = effect.B117
    )
  )
}

#' Fit brms model for B.1.1.7 analysis
#' 
#' @param model brm model formula
#' @param data data to be used
#' @param fn names of the model
#' @param adapt_delta parameter for estimation in Stan
#' @return brms fit object 
fit_B117_model = function(model,data,fn,adapt_delta = .8) {
  fn = here(paste0("B117fits/",fn,".Rdata"))
  if (file.exists(fn)) {
    load(fn)
  } else {
    Bfit = 
      brm(model,
          data = data,
          backend = "cmdstanr",
          iter = 3500,
          warmup = 1000,
          control = list(adapt_delta = adapt_delta))
    
    draws = as_draws(Bfit$fit)
    sampler_params = 
      nuts_params(Bfit) %>% 
      data.table() %>% 
      dcast(Chain + Iteration ~ Parameter, value.var = "Value")
    save(Bfit,draws,sampler_params, file = fn)
  }
  return(Bfit)
}

#' Initial values for estimation of time course model in Stan
#' 
#' @param datalist list with data for Stan model
#' @return List with initial values
make_TC_inits = 
  function(datalist) {
    list(beta_sweight_mu = rnorm(1,8,.5),
         log_slope_up_mu = rnorm(1,.5,.025),
         log_slope_down_mu = rnorm(1,-1.75,.05),
         intercept_sigma = runif(1,.0,.05),
         slope_up_sigma = runif(1,.0,.05),
         slope_down_sigma = runif(1,.0,.05),
         intercept_raw = rnorm(datalist$G,0,.5),
         slope_down_raw = rnorm(datalist$G,0,.5),
         slope_up_raw = rnorm(datalist$G,0,.5),
         centre_sigma = runif(1,.1,.6),
         centre_raw = rnorm(datalist$N_centres,0,.5),
         intercept_PCR = rnorm(1,0,.5),
         int_centre1_sigma = runif(1,.1,.5),
         int_centre1_raw = rnorm(datalist$N_centre1,0,.5),
         b_shift = rep(0,datalist$G),
         # b_shift = do.call(c,lapply(1:datalist$G, 
         #                             function(x) 
         #                               ifelse(x %in% datalist$T1_neg_idx,-5,5) + runif(1,-1,1))),
         shift_centre1_sigma = runif(1,.1,.5), #
         shift_centre1_raw = array(rnorm(datalist$N_centre1,0,.5),dim = 1),
         a_raw_intercept_Age = rnorm(datalist$num_basis_Age,0,.05),
         tau_intercept_Age = runif(1,0,.05),
         a_raw_slope_up_Age = rnorm(datalist$num_basis_Age,0,.05),
         tau_slope_up_Age = runif(1,0,.05),
         betaPGH_slope_up = rnorm(ncol(datalist$X_PGH),0,.05),
         a_raw_slope_down_Age = rnorm(datalist$num_basis_Age,0,.05),
         tau_slope_down_Age = runif(1,0,.05),
         betaPGH_slope_down = rnorm(ncol(datalist$X_PGH),0,.05),
         a_intercept_Age = rnorm(0,.5), #
         a0_intercept_Age = runif(1,0,.1),
         log_intercept_mu = rnorm(1,log(7),.5),
         betaPGH_intercept = rnorm(ncol(datalist$X_PGH),0,.01),
         a_slope_down_Age = rnorm(0,.5),
         slope_down_ld_centre_sigma = runif(1,.01,.05),
         slope_down_ld_centre_raw = rnorm(datalist$N_ld_centre,0,.5),
         sigma_mu = rnorm(1,0,.5),
         sigma_sigma = runif(1,.1,.2),
         sigma_raw = rnorm(datalist$G,0,.5),
         imp_neg = array(runif(datalist$N_NegTests,-1,datalist$imputation_limit),
                         dim = c(1,datalist$N_NegTests)),
         a0_slope_up_Age = rnorm(1,0,.05),
         a0_slope_down_Age = rnorm(1,0,.05)
    )
  }

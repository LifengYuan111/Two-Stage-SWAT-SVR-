# RCW two-stage SVR: five-fold expanding-window temporal CV
# R >= 4.2; packages dplyr, lubridate, e1071, readr, FlowScreen
# Stage 1: Eckhardt baseflow proxy; Stage 2: SWAT+ flow + BF proxy
# Tuning uses calibration only; validation is held out.
suppressPackageStartupMessages({
 library(dplyr); library(lubridate); library(e1071)
 library(readr); library(FlowScreen)
})
set.seed(123)

root <- "C:/Users/lyccf/University of Missouri/Fort Leonard Wood - Ogrp - General/Lifeng_files/Manuscript/ENVSOFT-D-26-04122/The revised version"
monthly_file <- file.path(root,"RCW_Monthly_Matched_Three_Models_With_Baseflow.csv")
daily_file <- file.path(root,"RCW_Daily_Matched_Three_Models_With_Baseflow.csv")
raw_file <- "C:/Users/lyccf/University of Missouri/Fort Leonard Wood - Ogrp - General/Lifeng_files/FLW calibration/USGS_Hydrologic_WQ_Data.csv"
out <- file.path(root,"SVR_TwoStage_Optimization_5Fold")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(f in c(monthly_file,daily_file,raw_file))
 if(!file.exists(f)) stop("Missing input: ",f)
area <- c("06928300"=427.6,"06928320"=25.2,"06928330"=13.4,
          "06928380"=21.7,"06928420"=706.6)
canon <- function(x) sprintf("%08d",as.integer(gsub("[^0-9]","",as.character(x))))
raw <- read.csv(raw_file,check.names=FALSE,colClasses=c(site_no="character")) %>%
 mutate(site_no=canon(site_no),Date=as.Date(Date),
        q=as.numeric(X_00060_00003)*0.0283168,
        flag=as.character(X_00060_00003_cd)) %>%
 filter(as.character(parameter_cd)%in%c("00060","60"),
        is.finite(q),q>=0,grepl("A",flag)) %>%
 select(site_no,Date,q) %>% distinct(site_no,Date,.keep_all=TRUE)

# Non-overlapping forward test windows; training is strictly earlier.
# d must contain CALIBRATION rows only.
folds <- function(d,k=5,initial=.40) {
 dates<-sort(unique(as.Date(d$Date))); n<-length(dates)
 if(n<k+10) stop("Too few distinct dates for ",k," folds")
 cuts<-unique(floor(seq(initial*n,n,length.out=k+1)))
 if(length(cuts)!=k+1 || any(diff(cuts)<2)) stop("Invalid fold boundaries")
 lapply(seq_len(k),function(f){
   i<-cuts[f]; j<-cuts[f+1]
   tr<-which(d$Date<=dates[i])
   te<-which(d$Date>dates[i] & d$Date<=dates[j])
   if(!length(tr)||!length(te)) stop("Empty fold")
   list(train=tr,test=te,train_end=dates[i],test_end=dates[j])
 })
}
# Compute Eckhardt target only up to a cutoff; zero flows are retained.
get_bf <- function(d,cutoff){
 d<-d %>% filter(Date<=cutoff) %>% arrange(Date)
 if(nrow(d)<3) return(tibble(Date=as.Date(character()),bf=numeric()))
 grid<-tibble(Date=seq(min(d$Date),max(d$Date),by="day")) %>% left_join(d,by="Date")
 v<-grid$q; ix<-which(is.finite(v))
 if(length(ix)<3) return(tibble(Date=grid$Date,bf=rep(NA_real_,length(v))))
 filled<-approx(ix,v[ix],xout=seq_along(v),rule=2)$y
 b<-as.numeric(FlowScreen::bf_eckhardt(discharge=filled,a=.98,BFI=.80))
 tibble(Date=grid$Date,bf=ifelse(is.finite(v),pmax(0,pmin(b,v)),NA_real_))
}
baseflow_targets <- function(cutoff,scale){
 if(scale=="monthly") cutoff<-ceiling_date(as.Date(cutoff),"month")-days(1)
 z<-bind_rows(lapply(split(raw,raw$site_no),function(d){
   b<-get_bf(d,cutoff); b$site_no<-d$site_no[1]; b
 }))
 if(scale=="monthly")
   z<-z %>% mutate(Date=floor_date(Date,"month")) %>%
    group_by(site_no,Date) %>%
    summarise(bf=if(sum(is.finite(bf))>=10) mean(bf,na.rm=TRUE)
              else NA_real_,.groups="drop")
 z
}
metrics <- function(y,p){
 ok<-is.finite(y)&is.finite(p);y<-y[ok];p<-p[ok]
 if(length(y)<3 || sum((y-mean(y))^2)<1e-12)
   return(c(n=length(y),NSE=NA_real_,RMSE=NA_real_,PBIAS=NA_real_))
 c(n=length(y),NSE=1-sum((y-p)^2)/sum((y-mean(y))^2),
   RMSE=sqrt(mean((y-p)^2)),
   PBIAS=if(sum(y)>0)100*sum(p-y)/sum(y) else NA_real_)
}
features <- function(d,scale){
 # WARNING: lag() on sparse matched dates can jump over gaps.
 # Only allow a lag when dates are exactly consecutive at the scale.
 d<-d %>% arrange(site_no,Date) %>% group_by(site_no) %>%
 mutate(month=month(Date),sin12=sin(2*pi*month/12),
        cos12=cos(2*pi*month/12)) %>% ungroup()
 # Explicitly nullify lags when calendar periods are not adjacent.
 if(scale=="monthly"){
   d<-d %>% group_by(site_no) %>%
    mutate(sim1=if_else((year(Date)*12+month(Date))-
                (year(lag(Date))*12+month(lag(Date)))==1,lag(sim_q),NA_real_),
           sim2=if_else((year(Date)*12+month(Date))-
                (year(lag(Date,2))*12+month(lag(Date,2)))==2,lag(sim_q,2),NA_real_),
           p1=if_else((year(Date)*12+month(Date))-
                (year(lag(Date))*12+month(lag(Date)))==1,lag(precip_mm),NA_real_),
           p2=if_else((year(Date)*12+month(Date))-
                (year(lag(Date,2))*12+month(lag(Date,2)))==2,lag(precip_mm,2),NA_real_)) %>% ungroup()
 } else {
   d<-d %>% group_by(site_no) %>%
    mutate(sim1=if_else(as.integer(Date-lag(Date))==1,lag(sim_q),NA_real_),
           sim2=if_else(as.integer(Date-lag(Date,2))==2,lag(sim_q,2),NA_real_),
           p1=if_else(as.integer(Date-lag(Date))==1,lag(precip_mm),NA_real_),
           p2=if_else(as.integer(Date-lag(Date,2))==2,lag(precip_mm,2),NA_real_)) %>% ungroup()
 }
 d$p3<-d$precip_mm+d$p1+d$p2
 d$area_km2<-unname(area[d$site_no])
 for(v in c("sim_q","precip_mm","sim1","sim2","p1","p2","p3","area_km2"))
   d[[paste0("x_",v)]]<-log1p(pmax(d[[v]],0))
 bfcols<-paste0("x_",c("sim_q","precip_mm","sim1","sim2",
                         "p1","p2","p3","area_km2"))
 if("api"%in%names(d)){d$x_api<-log1p(pmax(d$api,0));bfcols<-c(bfcols,"x_api")}
 bfcols<-c(bfcols,"sin12","cos12")
 list(data=d,bfcols=bfcols)
}
fit_predict <- function(train,test,cols,target,params){
 ok<-complete.cases(train[,c(cols,target),drop=FALSE]);train<-train[ok,,drop=FALSE]
 pred<-rep(NA_real_,nrow(test))
 if(nrow(train)<15)return(pred)
 mu<-vapply(train[,cols,drop=FALSE],mean,numeric(1))
 sdv<-vapply(train[,cols,drop=FALSE],sd,numeric(1))
 sdv[!is.finite(sdv)|sdv==0]<-1
 xt<-sweep(sweep(as.matrix(train[,cols,drop=FALSE]),2,mu),2,sdv,"/")
 xv<-sweep(sweep(as.matrix(test[,cols,drop=FALSE]),2,mu),2,sdv,"/")
 valid<-complete.cases(xv);if(!any(valid))return(pred)
 model<-e1071::svm(x=xt,y=log1p(pmax(train[[target]],0)),
    type="eps-regression",kernel="radial",epsilon=params$epsilon,
    gamma=params$gamma,cost=params$cost,scale=FALSE)
 pred[valid]<-pmax(0,expm1(as.numeric(predict(model,xv[valid,,drop=FALSE]))))
 pred
}
search_params <- function(base,n){
 bind_rows(as_tibble(base),
  tibble(epsilon=2^runif(n,-9,0),gamma=2^runif(n,-9,2),
         cost=2^runif(n,-2,9))) %>% distinct()
}
# Stage 2 training requires genuine OOF Stage 1 predictions.
oof_baseflow_scaled <- function(tr,bfcols,bestbf,scale,k=5){
 ff<-folds(tr,k=k); outbf<-rep(NA_real_,nrow(tr))
 for(g in ff){
   a<-tr[g$train,,drop=FALSE];z<-tr[g$test,,drop=FALSE]
   targets<-baseflow_targets(max(a$Date),scale)
   a<-a %>% select(-any_of("bf")) %>% left_join(targets,by=c("site_no","Date"))
   outbf[g$test]<-fit_predict(a,z,bfcols,"bf",bestbf)
 }
 outbf
}
run_scale <- function(scale, file, run_id = 1L, seed = 123L) {
 set.seed(seed)
 run_dir <- file.path(out,sprintf("Run_%02d",run_id))
 dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
 message("Starting ",scale," five-fold tuning")
 d<-read_csv(file,show_col_types=FALSE,
             col_types=cols(site_no=col_character(),.default=col_guess()))
 stopifnot(all(c("site_no","period","obs_q","sim_q","precip_mm","pred_svr")%in%names(d)))
 d$site_no<-canon(d$site_no)
 d$Date<-as.Date(if(scale=="monthly")d$Month else d$Date)
 if(anyDuplicated(d[c("site_no","Date")]))stop("Duplicate site/date: ",scale)
 if(!all(d$period%in%c("calibration","validation")))stop("Invalid period")
 if(any(!is.finite(d$obs_q))) stop("Nonfinite observed flow: ", scale)
 tmp<-features(d,scale);d<-tmp$data;bfcols<-tmp$bfcols
 cal<-d %>% filter(period=="calibration")
 val<-d %>% filter(period=="validation")
 fd<-folds(cal,k=5)
 bf_orig<-if(scale=="monthly")list(epsilon=.125,gamma=.0625,cost=16) else
                                   list(epsilon=.125,gamma=.25,cost=8)
 q_orig<-if(scale=="monthly")list(epsilon=.00390625,gamma=.5,cost=32) else
                                  list(epsilon=.25,gamma=.03125,cost=32)
 bfgrid<-search_params(bf_orig,if(scale=="monthly")48 else 24)
 # Precompute fold-specific targets, including evaluation target through fold end.
 bfcache<-lapply(fd,function(f){
   tr<-cal[f$train,,drop=FALSE];te<-cal[f$test,,drop=FALSE]
   tr<-tr %>% left_join(baseflow_targets(max(tr$Date),scale),by=c("site_no","Date"))
   te<-te %>% left_join(baseflow_targets(max(te$Date),scale),by=c("site_no","Date"))
   list(train=tr,test=te)
 })
 bf_scores<-bind_rows(lapply(seq_len(nrow(bfgrid)),function(i){
   p<-as.list(bfgrid[i,]);yy<-numeric();pp<-numeric()
   for(f in bfcache){
     pred<-fit_predict(f$train,f$test,bfcols,"bf",p)
     yy<-c(yy,f$test$bf);pp<-c(pp,pred)
   }
   tibble(bfgrid[i,],cv_NSE=unname(metrics(yy,pp)["NSE"]))
 })) %>% arrange(desc(cv_NSE))
 if(!is.finite(bf_scores$cv_NSE[1]))stop("No valid Stage 1 CV score: ",scale)
 bestbf<-as.list(bf_scores[1,c("epsilon","gamma","cost")])
 # For each outer fold, make inner five-fold OOF Stage 1 predictors.
 cache<-lapply(fd,function(f){
   tr<-cal[f$train,,drop=FALSE];te<-cal[f$test,,drop=FALSE]
   tr<-tr %>% left_join(baseflow_targets(max(tr$Date),scale),by=c("site_no","Date"))
   tr$bf_hat<-oof_baseflow_scaled(tr,bfcols,bestbf,scale,k=5)
   te$bf_hat<-fit_predict(tr,te,bfcols,"bf",bestbf)
   tr$x_bf_hat<-log1p(pmax(tr$bf_hat,0))
   te$x_bf_hat<-log1p(pmax(te$bf_hat,0))
   list(train=tr,test=te)
 })
 flowcols<-c("x_sim_q","x_bf_hat")
 qgrid<-search_params(q_orig,if(scale=="monthly")48 else 24)
 q_scores<-bind_rows(lapply(seq_len(nrow(qgrid)),function(i){
   p<-as.list(qgrid[i,]);yy<-numeric();pp<-numeric()
   for(f in cache){
     pred<-fit_predict(f$train,f$test,flowcols,"obs_q",p)
     yy<-c(yy,f$test$obs_q);pp<-c(pp,pred)
   }
   tibble(qgrid[i,],cv_NSE=unname(metrics(yy,pp)["NSE"]))
 })) %>% arrange(desc(cv_NSE))
 if(!is.finite(q_scores$cv_NSE[1]))stop("No valid Stage 2 CV score: ",scale)
 bestq<-as.list(q_scores[1,c("epsilon","gamma","cost")])
 cal<-cal %>% left_join(baseflow_targets(max(cal$Date),scale),by=c("site_no","Date"))
 cal$bf_hat<-oof_baseflow_scaled(cal,bfcols,bestbf,scale,k=5)
 val$bf_hat<-fit_predict(cal,val,bfcols,"bf",bestbf)
 cal$x_bf_hat<-log1p(pmax(cal$bf_hat,0))
 val$x_bf_hat<-log1p(pmax(val$bf_hat,0))
 val$pred_new<-fit_predict(cal,val,flowcols,"obs_q",bestq)
 old<-metrics(val$obs_q,val$pred_svr)
 new<-metrics(val$obs_q,val$pred_new)
 write_csv(bf_scores,file.path(run_dir,paste0(scale,"_stage1_search.csv")))
 write_csv(q_scores,file.path(run_dir,paste0(scale,"_stage2_search.csv")))
 write_csv(val %>% select(site_no,Date,period,obs_q,sim_q,pred_svr,bf_hat,pred_new),
           file.path(run_dir,paste0(scale,"_validation_predictions.csv")))
 write_csv(tibble(scale=scale,model=c("original_hybrid","new_hybrid"),
     NSE=c(old["NSE"],new["NSE"]),RMSE=c(old["RMSE"],new["RMSE"]),
     PBIAS=c(old["PBIAS"],new["PBIAS"])),
     file.path(run_dir,paste0(scale,"_validation_comparison.csv")))
 write_csv(tibble(stage=c("baseflow","streamflow"),
     epsilon=c(bestbf$epsilon,bestq$epsilon),
     gamma=c(bestbf$gamma,bestq$gamma),
     cost=c(bestbf$cost,bestq$cost)),
     file.path(run_dir,paste0(scale,"_selected_parameters.csv")))
 summary_row <- tibble(
   run_id=run_id,seed=seed,scale=scale,
   bf_epsilon=bestbf$epsilon,bf_gamma=bestbf$gamma,bf_C=bestbf$cost,
   bf_cv_NSE=bf_scores$cv_NSE[1],
   flow_epsilon=bestq$epsilon,flow_gamma=bestq$gamma,flow_C=bestq$cost,
   flow_cv_NSE=q_scores$cv_NSE[1],
   validation_NSE=unname(new["NSE"]),
   validation_RMSE=unname(new["RMSE"]),
   validation_PBIAS=unname(new["PBIAS"]),
   original_validation_NSE=unname(old["NSE"])
 )
 write_csv(summary_row,file.path(run_dir,paste0(scale,"_run_summary.csv")))
 message(scale," run ",run_id,": CV NSE=",round(summary_row$flow_cv_NSE,3),
         "; validation NSE=",round(summary_row$validation_NSE,3))
 return(summary_row)
}
# 20 reproducible randomized searches per scale. No validation-driven selection.
# Stop early if an outdated function definition was loaded.
stopifnot(all(c("scale", "file", "run_id", "seed") %in% names(formals(run_scale))))
n_runs <- as.integer(Sys.getenv("RCW_N_RUNS", "50"))
if (!is.finite(n_runs) || n_runs < 1L) stop("RCW_N_RUNS must be a positive integer")
all_results <- list()
for(run_id in seq_len(n_runs)) {
 for(scale in c("monthly","daily")) {
   input <- if(scale=="monthly") monthly_file else daily_file
   seed <- 123000L + 100L*run_id + if(scale=="daily") 1L else 0L
   result <- tryCatch(run_scale(scale,input,run_id=run_id,seed=seed),
     error=function(e){
       message("FAILED: ",scale," run ",run_id,": ",conditionMessage(e))
       tibble(run_id=run_id,seed=seed,scale=scale,
              error=conditionMessage(e))
     })
   all_results[[length(all_results)+1L]] <- result
   write_csv(bind_rows(all_results),file.path(out,"All_50_Runs_Progress.csv"))
 }
}
all_summary <- bind_rows(all_results)
write_csv(all_summary,file.path(out,"All_50_Runs_Summary.csv"))
valid <- if ("flow_cv_NSE" %in% names(all_summary)) {
  all_summary %>% filter(is.finite(flow_cv_NSE))
} else {
  tibble()
}
if(nrow(valid) > 0L) {
 ranking <- valid %>% group_by(scale) %>%
   arrange(desc(flow_cv_NSE),.by_group=TRUE) %>%
   mutate(cv_rank=row_number(),best_by_CV=cv_rank==1L) %>% ungroup()
 write_csv(ranking,file.path(out,"All_Runs_Ranked_By_CV_NSE.csv"))
 best <- ranking %>% filter(best_by_CV)
 write_csv(best,file.path(out,"BEST_Parameters_By_CV_NSE.csv"))
 print(best)
} else warning("No successful searches with finite CV NSE; check All_50_Runs_Summary.csv")
writeLines(capture.output(sessionInfo()),file.path(out,"R_session_info.txt"))
message("20-run search complete. Results saved to: ",out)

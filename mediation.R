assocMediLm = function(x,m,y,cov,dat, mediMethod='gformula',ciInfer='bootstrap', nboot=10){
  subDat = dplyr::select(dat, all_of(c(x,m,y,cov)))
  subDat = na.omit(subDat)
  mreg = purrr::map_chr(m, function(x){
    if(inherits(dat[[x]], c('numeric'))){
      mreg = 'linear'
    }else if(inherits(dat[[x]], c('factor'))){
      if(length(levels(dat[[x]]))==2){
        mreg = 'logistic'
      }else{
        stop('only have one level')
      }
    }
    mreg
  })
  yreg = purrr::map_chr(y, function(x){
    if(inherits(dat[[x]], c('numeric'))){
      mreg = 'linear'
    }else if(inherits(dat[[x]], c('factor'))){
      if(length(levels(dat[[x]]))==2){
        mreg = 'logistic'
      }else{
        stop('only have one level')
      }
    }
    mreg
  })
  est = do.call('cmest', list(data = subDat, model = mediMethod, outcome = y, exposure = x,
                              mediator = m, EMint = F, basec =  cov, 
                              mreg = as.list(mreg), yreg = yreg, astar = 0, a = 1, mval = as.list(rep(0,length(m))), 
                              inference = ciInfer, nboot = nboot))
}

tidyMediLm = function(fit){
  estSumm = summary(fit)
  res = estSumm$summarydf  %>% rownames_to_column('term')
  colnames(res) = c('term', 'estimate', 'se', 'lwr', 'upr', 'p')
  if(fit$call$yreg %in% c('logistic')){
    indexName = c("Rte", "Rpnde", "Rtnie", "pm")
  }else{
    indexName = c("te", "pnde", "tnie", "pm")
  }
  res = res[match(indexName, res$term),]
  res$term = c('te', 'de', 'ide', 'pm')
  res = tidyr::pivot_wider(res, names_from = 'term', values_from = c('estimate', 'se', 'lwr', 'upr', 'p'))
  res$x = estSumm$call$exposure
  res$m = paste0(estSumm$call$mediator, collapse = ', ')
  res$y = estSumm$call$outcome
  res$yreg = estSumm$call$yreg
  res
}

applymediLm = function(xl, ml, yl, covl, dat, mediMethod='gformula',ciInfer='bootstrap', nboot=10, expand = F, onlySummary = F, cores=4, outdir=NULL){
  varlist = expandVars(xl, ml, yl, covl, colname = c('x', 'm', 'y', 'cov'), expand = expand)
  progressr::with_progress({
    future::plan(future::multisession, workers = cores)
    pb = progressr::progressor(steps  = 100)
    res = furrr::future_map(1:nrow(varlist), function(i, xl, ml, yl, covl, varlist, dat, mediMethod, ciInfer, nboot, onlySummary, pb){
      library(CMAverse)
      x=unlist(varlist[i,]$x)
      m=unlist(varlist[i,]$m)
      y=unlist(varlist[i,]$y)
      cov=unlist(varlist[i,]$cov)
      if(i %% round(nrow(varlist)/100) == 0 | i == nrow(varlist)) pb()
      if(!is.null(outdir)){
        if(!dir.exists(outdir)){dir.create(outdir)}
        fp = paste0(digest::digest(paste(x,m,y, sep='_')),'.RDS')
        fp = file.path(outdir, fp)
        if(file.exists(fp)){res = readRDS(fp);return(res)}
        
      }
      fit = assocMediLm(x=x,m=m,y=y, cov=cov, dat=dat, mediMethod=mediMethod, ciInfer=ciInfer, nboot=nboot)
      tidyRes =tidyMediLm(fit)
      res = list('Summary'=tidyRes, x=x, m=paste0(m, collapse = ', '), y=y, cov=list(cov))
      if(!onlySummary){res$fit = fit}
      if(!is.null(outdir)){
        write_rds(res, fp)
      }
      res
    }, xl = xl, ml=ml, yl = yl, covl = covl, varlist = varlist, dat = dat, mediMethod=mediMethod, ciInfer=ciInfer, nboot = nboot, onlySummary = onlySummary, pb = pb)
    future::plan(future::sequential)
  })
  res
}

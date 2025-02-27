assocttest = function(x, y, dat, paired = F){
  testRes  = t.test(myfml(x = x, cov = NULL, y = y, rnd=NULL), paired = paired, data = dat)
  testRes
}
tidyttest = function(fit, x, y){
  res = broom::tidy(fit)
  res$x = x
  res$y = y
  res
}
applyttest = function(xl, yl, dat, expand = F, paired = F, cores, onlySummary = T, outdir=NULL){
  varlist = expandVars(xl, yl, covl= NULL, colname = c('x', 'y', 'cov'), expand=expand)
  handlers(handler_progress(interval = 5))
  progressr::with_progress({
    future::plan(future::multisession, workers = cores)
    pb = progressr::progressor(steps = nrow(varlist))
    res = furrr::future_map(1:nrow(varlist), function(i, varlist, dat, paired, pb){
      library(broom)
      x=unlist(varlist[i,]$x)
      y=unlist(varlist[i,]$y)
      fulln = paste0('temp_', paste0(c(y,x), collapse = '_'), '_CVR2.rds')
      if(is.null(outdir)){
        dat = na.omit(dplyr::select(dat, all_of(c(x,y))))
        fit = assocttest(x, y, dat, paired)
        tidyRes = tidyttest(fit, x, y)
        res = list('Summ'=tidyRes, x=x, y=y)
        if(!onlySummary){res$fit = fit}
      }else{
        files = list.files(outdir, include.dirs = F)
        hasDone = any(str_detect(files, fulln))
        if(hasDone){
          res = readRDS(file.path(outdir, fulln))
        }else{
          dat = na.omit(dplyr::select(dat, all_of(c(x,y))))
          fit = assocttest(x, y, dat, paired)
          tidyRes = tidyttest(fit, x, y)
          res = list('Summ'=tidyRes, x=x, y=y)
          if(!onlySummary){res$fit = fit}
          saveRDS(res, file.path(outdir, fulln))
        }
      }
      pb()
      res
    }, varlist = varlist, dat = dat, paired = paired, pb = pb)
    future::plan(future::sequential)
  })
  res
}

myfml = function(x, cov, y, rnd, stringAsFormula=T){
  x[is.na(x)] = NULL
  y[is.na(y)] = NULL
  cov[is.na(cov)] = NULL
  rnd[is.na(rnd)] = NULL
  emptystr2NULL = function(x){if(x==''|is.na(x)) NULL else x}
  ystr = emptystr2NULL(paste0(y, collapse = ' + '))
  fusex = na.omit(c(x, cov))
  xstr = emptystr2NULL(paste0(fusex, collapse = ' + '))
  rndstr = emptystr2NULL(paste0(rnd, collapse = ' + '))
  right = emptystr2NULL(paste0(c(xstr, rndstr), collapse = ' + '))
  exprstr = paste0(ystr, ' ~ ', right, collapse = '')
  if (stringAsFormula) as.formula(exprstr) else exprstr
}

assocLm = function(x, y, cov, dat, family){
  fml = myfml(x=x, cov=cov, y=y, rnd=NULL)
  fit = do.call('glm', args = list(formula = fml, family = family, data=dat))
  # fit = glm(formula = fml, family = family, data=dat)
  fit
}

tidyLm = function(fit, x=NULL, y=NULL, calR2=F){
  # T test
  tres = broom::tidy(fit)
  tres$nobservation = length(fit$residuals)
  tres$x = x
  tres$y = y
  # F test
  anovares = car::Anova(fit, type=2)
  anovares = rownames_to_column(as.data.frame(anovares), 'term')
  # print(colnames(anovares))
  anovares = dplyr::rename(anovares, `p.value`=`Pr(>Chisq)`)
  anovares$nobservation = length(fit$residuals)
  anovares$x = x
  anovares$y = y
  if(calR2){
    rfit = update(fit, as.formula(paste0('~. -', x)))
    r2 = rsq::rsq.partial(fit, rfit)$partial.rsq
    tres$R2 = r2
    anovares$R2 = r2
  }
  ret = list('T'=tres, 'F'=anovares)
  ret
}

expandVars = function(..., colname = NULL, expand=F){
  varlist = list(...)
  if(expand){
    varlist = expand.grid(purrr::map(varlist, function(x){if(is_empty(x)){as.list(c(NA))}else{x}}))
    varlist = rlist::list.cbind(purrr::map(varlist, function(x){purrr::map(x, function(y){if(all(is.na(y))){NULL}else{y}})}))
  }else{
    varlist = as.data.frame(rlist::list.cbind(varlist))
  }
  colnames(varlist) = if(all(is.null(colname))){paste0('var', 1:ncol(varlist))} else {colname}
  varlist
}

applyLm = function(xl, yl, covl, dat, family='gaussian', expand = F, onlySummary = T, calR2=F, cores=4, outdir=NULL){
  if (!is.null(outdir) && !dir.exists(outdir)) dir.create(outdir)
  varlist = expandVars(xl, yl, covl, colname = c('x', 'y', 'cov'), expand=expand)
  handlers(handler_progress(interval = 5))
  progressr::with_progress({
    future::plan(future::multisession, workers = cores)
    pb = progressr::progressor(steps = nrow(varlist))
    res = furrr::future_map(1:nrow(varlist), function(i, varlist, dat, family, onlySummary, calR2, pb){
      library(broom)
      x=unlist(varlist[i,]$x)
      y=unlist(varlist[i,]$y)
      cov=unlist(varlist[i,]$cov)
      fulln = paste0('temp_', paste0(c(y,x), collapse = '_'), '_CVR2.rds')
      if(is.null(outdir)){
        dat = na.omit(dplyr::select(dat, all_of(c(x,y,cov))))
        fit = assocLm(x, y, cov, dat, family)
        tidyRes = tidyLm(fit, x=x, y=y, calR2 = calR2)
        res = list('T'=tidyRes$`T`, 'F'=tidyRes$`F`, x=x, y=y, cov=list(cov))
        if(!onlySummary){res$fit = fit}
      }else{
        files = list.files(outdir, include.dirs = F)
        hasDone = any(str_detect(files, fulln))
        if(hasDone){
          res = readRDS(file.path(outdir, fulln))
        }else{
          dat = na.omit(dplyr::select(dat, all_of(c(x,y,cov))))
          fit = assocLm(x, y, cov, dat, family)
          tidyRes = tidyLm(fit, x=x, y=y, calR2 = calR2)
          res = list('T'=tidyRes$`T`, 'F'=tidyRes$`F`, x=x, y=y, cov=list(cov))
          if(!onlySummary){res$fit = fit}
          saveRDS(res, file.path(outdir, fulln))
        }
      }
      pb()
      res
    }, varlist = varlist, dat = dat, family = family, onlySummary = onlySummary, calR2=calR2, pb = pb)
    future::plan(future::sequential)
  })
  res
}

assocLmm = function(x, y, cov, rnd, dat, family){
  require(lmerTest, quietly =T, warn.conflicts=F)
  fml = myfml(x=x, cov=cov, y=y, rnd=rnd)
  if (family == 'gaussian'){
    # fit =lmer(fml, data = dat, control = lmerControl(optimizer = 'Nelder_Mead', check.conv.singular = 'ignore'))
    fit = do.call('lmer', args = list(formula=fml, data = dat, control =  lmerControl(optimizer = 'Nelder_Mead', check.conv.singular = 'ignore')))
  }else{
    # fit =glmer(fml, family = family, data = dat, control = glmerControl(optimizer = 'Nelder_Mead', check.conv.singular = 'ignore'))
    fit = do.call('glmer', args = list(formula=fml, family = family, data = dat, control =  glmerControl(optimizer = 'Nelder_Mead', check.conv.singular = 'ignore')))
  }
  fit
}

tidyLmm = function(fit, x=NULL, y=NULL, calR2=F){
  library(broom)
  library(broom.mixed)
  library(lmerTest)
  # require(broom, quietly =T, warn.conflicts=F)
  # require(broom.mixed, quietly =T, warn.conflicts=F)
  # require(lmerTest, quietly =T, warn.conflicts=F)
  tres = broom::tidy(fit)
  summ = summary(fit)
  tres$nparticipant = summ$ngrps[!str_detect(names(summ$ngrps), ':')]
  tres$nobservation = nrow(fit@frame)
  tres$x = x
  tres$y = y
  anovares = car::Anova(fit, type=3)
  anovares = rownames_to_column(as.data.frame(anovares), 'term')
  anovares = dplyr::rename(anovares, p.value=`Pr(>Chisq)`)
  anovares$nparticipant = summ$ngrps[!str_detect(names(summ$ngrps), ':')]
  anovares$nobservation = length(summ$residuals)
  anovares$x = x
  anovares$y = y
  if(calR2){
    partR2res = partR2::partR2(fit, partvars = x, max_level = 1)
    tres$R2 = unlist(partR2res$R2[2,2])
    anovares$R2 = unlist(partR2res$R2[2,2])
  }
  ret = list('T'=tres, 'F'=anovares)
  ret
}

applyLmm = function(xl, yl, covl, rndl, dat, family='gaussian', expand = F, onlySummary = F, calR2=F, cores=4, outdir=NULL){
  if (!is.null(outdir) && !dir.exists(outdir)) dir.create(outdir)
  varlist = expandVars(xl, yl, covl, rndl, colname = c('x', 'y', 'cov', 'rnd'), expand = expand)
  handlers(handler_progress(interval = 5))
  progressr::with_progress({
    future::plan(future::multisession, workers = cores)
    pb = progressr::progressor(steps = nrow(varlist))
    # res = furrr::future_map(1:nrow(varlist), function(i, xl, yl, covl, rndl, varlist, dat, family, onlySummary, calR2, pb){
    res = furrr::future_map(1:nrow(varlist), function(i, varlist, dat, family, onlySummary, calR2, outdir, pb){
      require(broom, quietly =T, warn.conflicts=F)
      require(lmerTest, quietly =T, warn.conflicts=F)
      x=unlist(varlist[i,]$x)
      y=unlist(varlist[i,]$y)
      cov=unlist(varlist[i,]$cov)
      rnd = unlist(varlist[i,]$rnd)
      fulln = paste0('temp_', paste0(c(y,x), collapse = '_'), '_CVR2.rds')
      if(is.null(outdir)){
        fit = assocLmm(x=x, y=y, cov=cov, rnd=rnd, dat=dat, family=family)
        tidyRes = tidyLmm(fit, x=x, y=y, calR2=calR2)
        res = list('T'=tidyRes$`T`, 'F'=tidyRes$`F`, x=x, y=y, cov=list(cov), rnd=list(rnd))
        if(!onlySummary){res$fit = fit}
      }else{
        files = list.files(outdir, include.dirs = F)
        hasDone = any(str_detect(files, fulln))
        if(hasDone){
          res = readRDS(file.path(outdir, fulln))
        }else{
          fit = assocLmm(x=x, y=y, cov=cov, rnd=rnd, dat=dat, family=family)
          tidyRes = tidyLmm(fit, x=x, y=y, calR2=calR2)
          res = list('T'=tidyRes$`T`, 'F'=tidyRes$`F`, x=x, y=y, cov=list(cov), rnd=list(rnd))
          if(!onlySummary){res$fit = fit}
          saveRDS(res, file.path(outdir, fulln))
        }
      }
      pb()
      res
      # }, xl = xl, yl = yl, covl = covl, rndl = rndl, varlist = varlist, dat = dat, family = family, onlySummary = onlySummary, calR2=calR2, pb = pb)
    }, varlist = varlist, dat = dat, family = family, onlySummary = onlySummary, calR2=calR2, outdir = outdir, pb = pb)
    future::plan(future::sequential)
  })
  res
}


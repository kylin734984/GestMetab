predictgbpCVR2 = function(x, y, dat, nfold, disable_neg = F, bs_num =10, seed=123, cores = 4, obj_id='pid', likelihood = "gaussian", search_params = list(), other_params = list(), num_try_random = 10, estimate_shap = T, nrounds=100, saveModel=T, bootParams = list('learning_rate'=0.01, 'max_depth'=5, 'feature_fraction'=0.8, 'num_leaves'=25, 'min_data_in_leaf'=15, 'n_estimators'=200, 'bagging_fraction'=0.9, 'metric'='mse')){
  gpbFitPredict = function(x, y, foldi, dobs, saveModel, saveData, dat, seed, cores, obj_id, likelihood, search_params, other_params, num_try_random, metric, nrounds){
    set.seed(seed)
    predVal = numeric(length = nrow(dat))
    fits = list()
    truey = dat[[y]]
    for (i in sort(unique(foldi))){
      print(paste0('fold ', i))
      trDat = dat[foldi != i, ]
      teDat = dat[foldi == i, ]
      trDat = if (nrow(trDat)==0) teDat else  trDat
      trDat = if (dobs)  trDat[trDat[[obj_id]] %in% sample(unique(trDat[[obj_id]]), replace = T),] else trDat
      trX = as.matrix(dplyr::select(trDat, all_of(x)))
      try = trDat[[y]]
      trGroup = dplyr::select(trDat, all_of(obj_id))
      teX = as.matrix(dplyr::select(teDat, all_of(x)))
      tey = teDat[[y]]
      teGroup = dplyr::select(teDat, all_of(obj_id))
      # trGpmodel <- GPModel(group_data = trGroup, likelihood = likelihood)
      trGpmodel <- GPModel(group_data = trGroup, likelihood = likelihood, gp_approx = 'vecchia', cov_function = 'wendland')
      trGpmodel$set_optim_params(params=list(optimizer_cov="nelder_mead", delta_rel_conv=1e-3)) # 针对大数据加速
      trdataset = gpb.Dataset(data = trX, label = try)
      if (all(purrr::map_dbl(search_params, length)==1)){
        print('test')
        fit <- gpb.train(data = trdataset, gp_model = trGpmodel, nrounds = nrounds, params = c(search_params, other_params), verbose = 0, 
                         objective = 'regression_l2', train_gp_model_cov_pars = F, num_threads = cores)
      }else{
        opt_params <- gpb.grid.search.tune.parameters(param_grid = search_params, params = other_params,
                                                      num_try_random = num_try_random, nfold = 3,
                                                      data = trdataset, gp_model = trGpmodel,
                                                      use_gp_model_for_validation = F, verbose_eval = 0,
                                                      nrounds = nrounds, early_stopping_rounds = 20, 
                                                      objective = 'regression_l2', num_threads = cores)
        fit <- gpb.train(data = trdataset, gp_model = trGpmodel, nrounds = opt_params$best_iter,
                         params = c(opt_params$best_params, other_params), verbose = 0, objective = 'regression_l2', 
                         train_gp_model_cov_pars = F, num_threads = cores)
      }
      pred_resp <- predict(fit, data = teX, group_data_pred = teGroup, predict_var = FALSE, pred_latent = FALSE, predict_disable_shape_check=TRUE)
      predVal[foldi == i] = pred_resp$response_mean
      fits = c(fits, list(fit))
    }
    res = list()
    res$truey = truey
    res$prediction = predVal
    res$fits = if(saveModel) fits else NULL
    res$data = if(saveData) dat else NULL
    res$R2 = mlr3measures::rsq(res$truey, res$prediction)
    res
  }
  foldi = foldIndexing(dat=dat, obj_id=obj_id, nfold=nfold, seed=seed)
  res = gpbFitPredict(x=x, y=y, foldi=foldi, dobs=F, saveModel=T, saveData=F, dat=dat, seed=seed, obj_id = obj_id, likelihood=likelihood, search_params=search_params, other_params=other_params, num_try_random=num_try_random, metric=metric, cores = cores, nrounds=nrounds)
  res$foldi = foldi
  res$boostrap_predictions = NULL
  res$boostrap_R2 = NULL
  res$shap_value = NULL
  if (estimate_shap){
    shapval = as.data.frame(matrix(NA, ncol = length(x), nrow = nrow(dat)))
    colnames(shapval) = x
    # shapval = as_tibble(matrix(numeric(nrow(dat)*ncol(dat)), ncol = ncol(dat)))
    for (i in unique(res$foldi)){
      subtraindat = as.matrix(dat[res$foldi == i, x])
      colnames(subtraindat) = x
      shapval[res$foldi == i, ] = SHAPforxgboost::shap.values(res$fits[[i]], subtraindat)$shap_score
    }
    res$shap_value = shapval
  }
  if ((disable_neg & res$R2<0) | bs_num<1) return(res)
  future::plan(future::multisession, workers = cores)
  progressr::with_progress({
    bspreds = furrr::future_map(1:bs_num, function(bsi, x, y, dat, foldi, seed, obj_id, likelihood, num_try_random){
      if(is.null(bootParams)){
        res = gpbFitPredict(x=x, y=y, foldi=foldi, dobs=T, saveModel=F, saveData=F, dat=dat, seed=seed+bsi, obj_id = obj_id, likelihood=likelihood, search_params=search_params, other_params=other_params, num_try_random=num_try_random, cores=1, nrounds=nrounds)
      }else{
        res = gpbFitPredict(x=x, y=y, foldi=foldi, dobs=T, saveModel=F, saveData=F, dat=dat, seed=seed+bsi, obj_id = obj_id, likelihood=likelihood, search_params=bootParams, other_params=other_params, num_try_random=num_try_random, cores=1, nrounds=nrounds)
      }
      
    }, x=x, y=y, dat=dat, foldi=foldi, seed=seed, obj_id = obj_id, likelihood=likelihood, num_try_random=num_try_random, .options = furrr_options(seed=T))
  })
  future::plan(future::sequential)
  res$boostrap_predictions = rlist::list.cbind(rlist::list.map(bspreds, prediction))
  res$boostrap_R2 = unlist(rlist::list.map(bspreds, R2))
  res$fits = if(saveModel) fits else NULL
  res
}

applyPredictGpbCVR2 = function(xl, yl, expand = F, dat, nfold, disable_neg=F, bs_num=10, seed=123, cores=4, obj_id='pid', likelihood = "gaussian", search_params = list(), other_params = list(), num_try_random = 10, estimate_shap=T, nrounds=100, saveModel=F, bootParams = NULL, outdir = NULL){
  if (!is.null(outdir) && !dir.exists(outdir)) dir.create(outdir)
  varlist = expandVars(xl, yl, colname = c('x', 'y'), expand = expand)
  future::plan(future::multisession, workers=cores)
  progressr::with_progress({
    pb = progressr::progressor(along = xl)
    res = furrr::future_map(1:nrow(varlist), function(i, varlist, dat, nfold, disable_neg, bs_num, seed, obj_id, likelihood, search_params, other_params, num_try_random, metric, estimate_shap, nrounds, outdir, saveModel, bootParams, pb){
      x=unlist(varlist[i,]$x)
      y=unlist(varlist[i,]$y)
      fulln = paste0(digest::digest(paste(x, y, sep='_')),'.RDS') #
      if(is.null(outdir)){
        res = predictgbpCVR2(x=x, y=y, dat=dat, nfold=nfold, disable_neg=disable_neg, bs_num=bs_num, seed=seed, cores=1, obj_id=obj_id, likelihood=likelihood, search_params=search_params, other_params=other_params, num_try_random=num_try_random, estimate_shap=estimate_shap, nrounds=nrounds, saveModel = saveModel, bootParams=bootParams)
        res$x = x
        res$y = y
      }else{
        files = list.files(outdir, include.dirs = F)
        hasDone = any(str_detect(files, fulln))
        if(hasDone){
          res = readRDS(file.path(outdir, fulln))
        }else{
          res = predictgbpCVR2(x=x, y=y, dat=dat, nfold=nfold, disable_neg=disable_neg, bs_num=bs_num, seed=seed, cores=1, obj_id=obj_id, likelihood=likelihood, search_params=search_params, other_params=other_params, num_try_random=num_try_random, estimate_shap=estimate_shap, nrounds=nrounds, saveModel = saveModel, bootParams=bootParams)
          res$x = x
          res$y = y
          saveRDS(res, file.path(outdir, fulln))
        }
      }
      pb()
      res
    }, varlist = varlist, dat=dat, nfold=nfold, disable_neg=disable_neg, bs_num=bs_num, seed=seed, obj_id=obj_id, likelihood=likelihood, search_params=search_params, other_params=other_params, num_try_random=num_try_random, estimate_shap=estimate_shap, nrounds=nrounds, saveModel=saveModel, bootParams=bootParams, outdir = outdir, pb = pb, .options = furrr_options(seed=T))
  })
  future::plan(future::sequential)
  if(!is.null(outdir)){
    files = list.files(outdir, include.dirs = F)
    fullns = map_chr(1:nrow(varlist), function(i, varlist){
      x=unlist(varlist[i,]$x)
      y=unlist(varlist[i,]$y)
      fulln = paste0(digest::digest(paste(x, y, sep='_')),'.RDS') # 因为x太多可能名字太长
    }, varlist=varlist)
    concatf = files[str_detect(files, paste0(fullns, collapse = '|'))]
    concatf = file.path(outdir,concatf)
    res = purrr::map(concatf, readRDS)
    saveRDS(res, file.path(outdir, 'allresults_CVR2.rds'))
  }
  res
}

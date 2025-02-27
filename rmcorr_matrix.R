cal_rep_corMat = function(individual_id, vars1, vars2=NULL, dat, cores=4, outdir=NULL){
  if (!is.null(outdir) && !dir.exists(outdir)) dir.create(outdir)
  dat[[individual_id]] = as.factor(dat[[individual_id]])
  if (!is.null(vars2)){
    ijs = as.data.frame(t(expand.grid(vars1, vars2)))
  }else{
    ijs = as.data.frame(combn(vars1, 2))
  }
  future::plan(future::multisession, workers = cores)
  progressr::with_progress({
    pb = progressr::progressor(steps = ncol(ijs))
    rhos = furrr::future_map_dfr(ijs, function(oneij, dat, individual_id, outdir, pb){
      fulln = paste0('temp_', paste0(c(oneij[1], oneij[2]), collapse = '_'), '_CVR2.rds')
      dat = dplyr::select(dat, all_of(c(oneij[1], oneij[2], individual_id)))
      if(is.null(outdir)){
        res = rmcorr::rmcorr(eval(parse(text = individual_id)), oneij[1], oneij[2], dat, CIs = 'analytic')
      }else{
        files = list.files(outdir, include.dirs = F)
        hasDone = any(str_detect(files, fulln))
        if(hasDone){
          res = readRDS(file.path(outdir, fulln))
        }else{
          res = rmcorr::rmcorr(eval(parse(text = individual_id)), oneij[1], oneij[2], dat, CIs = 'analytic')
          res$model <- NULL
          saveRDS(res, file.path(outdir, fulln))
        }
      }
      pb()
      data.frame('r' = res$r, 'df' = res$df, 'p' = res$p, 'ci' = paste(res$CI, sep='-', collapse = '-'), 
                 'var1' = res$vars[2],
                 'var2' = res$vars[3],
                 'cilevel' = res$CI.level)
    }, dat = dat, individual_id = individual_id, outdir = outdir, pb = pb)
  })
  future::plan(future::sequential)
  
  # obtain the correlation matrix
  ijs[3, ] = rhos$r
  ijs[4, ] = rhos$p
  if (is.null(vars2)){
    finalMatrix = matrix(nrow = length(vars1), ncol = length(vars1))
    dimnames(finalMatrix) = list(vars1, vars1)
    finalMatrixp = matrix(nrow = length(vars1), ncol = length(vars1))
    dimnames(finalMatrixp) = list(vars1, vars1)
    for (j in seq(1, ncol(ijs))){ finalMatrix[ijs[1,j], ijs[2,j]] = as.numeric(ijs[3,j])}
    for (oneX in vars1){finalMatrix[oneX, oneX] = 1}
    finalMatrix[lower.tri(finalMatrix)] = t(finalMatrix)[lower.tri(finalMatrix)]
    stopifnot(isSymmetric(finalMatrix))
    for (j in seq(1, ncol(ijs))){finalMatrixp[ijs[1,j], ijs[2,j]] = as.numeric(ijs[4,j])}
    for (oneX in vars1){finalMatrixp[oneX, oneX] = 0}
    finalMatrixp[lower.tri(finalMatrixp)] = t(finalMatrixp)[lower.tri(finalMatrixp)]
    stopifnot(isSymmetric(finalMatrixp))
  }else{
    finalMatrix = matrix(nrow = length(vars1), ncol = length(vars2))
    dimnames(finalMatrix) = list(vars1, vars2)
    for (j in seq(1, ncol(ijs))){finalMatrix[ijs[1,j], ijs[2,j]] = as.numeric(ijs[3,j])}
    finalMatrixp = matrix(nrow = length(vars1), ncol = length(vars2))
    dimnames(finalMatrixp) = list(vars1, vars2)
    for (j in seq(1, ncol(ijs))){finalMatrixp[ijs[1,j], ijs[2,j]] = as.numeric(ijs[4,j])}
  }
  res = list()
  res[['r']] = finalMatrix
  res[['p']] = finalMatrixp
  res[['summary']] = rhos
  return(res)
}
mymeta = function(estimate, se, idvars, datl, cores = 4){
  sub = function(estimate, se, idvars, dat, pb){
    suppressPackageStartupMessages(library(metafor)) 
    metaRes = do.call('rma', list(yi = dat[[estimate]], sei=dat[[se]], method = 'FE'))
    pb()
    metaSumm = tibble(beta_meta = as.numeric(metaRes$beta), se_meta = metaRes$se, beta_lci = metaRes$ci.lb, beta_uci = metaRes$ci.ub,
                      pval_meta = metaRes$pval, I2_meta = metaRes$I2, QEp = metaRes$QEp)
    metaSumm = dplyr::bind_cols(
      unique(dplyr::select(dat, all_of(idvars))),
      metaSumm)
  }
  future::plan(future::multisession, workers = cores)
  progressr::with_progress({
    pb = progressr::progressor(steps = length(datl))
    res = furrr::future_map_dfr(datl, sub, estimate=estimate, se =se, idvars =idvars, pb= pb, .options = furrr::furrr_options(seed = T))
  })
  future::plan(future::sequential)
  res
}
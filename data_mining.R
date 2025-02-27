imputeDat = function(dat, conImpute = min, catImpute = function(x){ifelse(is.na(x), '99', x)}){
  dat = dplyr::mutate_if(dat, is.numeric, function(x){as.numeric(Hmisc::impute(x, conImpute, na.rm = T))})
  dat = dplyr::mutate_if(dat, function(x){is.character(x)|is.factor(x)}, catImpute)
  return(dat)
}

namingout = function(prefixs, outsuffix, fmt=NULL){
  n1 = paste0(c(prefixs, outsuffix), collapse = '_')
  if (is.null(fmt)) n1 else paste(n1, fmt, sep='.')
}

SD_impu = function(x, sd_fold = 3){
  x[x<(mean(x, na.rm=T) - sd_fold*sd(x, na.rm = T))] = mean(x, na.rm=T) - sd_fold*sd(x, na.rm = T)
  x[x>(mean(x, na.rm=T) + sd_fold*sd(x, na.rm = T))] = mean(x, na.rm=T) + sd_fold*sd(x, na.rm = T)
  x
}
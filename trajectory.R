predictLmmTraj = function(x, y, pid, dat, d=2){
  # print('lmm')
  subdat = dplyr::select(dat, all_of(c(x,y,pid))) %>% na.omit()
  fits = purrr::map(1:d, function(i){
    require(splines, warn.conflicts=F)
    fml = myfml(x=paste0('bs(weeks, degree=', i, ')'), y = y, cov=NULL, rnd = paste0('(1|',pid, ')' ))
    lmerTest::lmer(fml, data = subdat)
  })
  opt = which.min(purrr::map_dbl(fits, AIC))
  fit = fits[[opt]]
  pred = tibble({{x}}:=4:40)
  pred[[y]] = as.numeric(predict(fit, newdata = pred, re.form=NA))
  pred
}



# trajectory plot ---------------------------------------------------------
p = ggplot(data = pdat, aes(x=x, y =y, color = stratafct)) +
  geom_line(stat = 'smooth', alpha = 1, linewidth = 0.5, se=F, method = 'gam') +
  labs(x = 'weeks', y = 'Z score', title = mapvalues(amet, metAnno$`COMP ID`, metAnno$compound, warn_missing = F)) + 
  scale_color_manual(values = c('0'='#1EB5BB', '1'='#EB746A'), breaks = c(0, 1), labels = c('Non HDP', 'HDP'), name = 'Prednisone') + 
  labs(color = stratafct) + 
  facet_wrap(~pop, scales = 'fixed', labeller = labeller(pop = popstripmap)) +
  theme_classic() +
  theme(legend.position = 'right', 
        text = element_text(size = 6), 
        strip.text = element_text(size = 6), 
        line = element_line(linewidth = 0.5),
        strip.background = element_rect(linewidth = 0.5))

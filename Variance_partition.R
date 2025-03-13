fit = lme4::lmer(value ~ income_class+residence+gestAge+para+edu_f + (1|pid), data = dat)
var_df = as.data.frame(insight::get_variance(fit))

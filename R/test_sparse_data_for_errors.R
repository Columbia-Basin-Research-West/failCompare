bad_time_dat=c(rep(5,3),rep(10,3))
my_e=quote(taglife.fn_weib3(tags.in = bad_time_dat,model.in = "weibull",tag.se = T))
my_e2=quote(taglife.fn_weib3(tags.in = bad_time_dat,model.in = "weibull",tag.se = F))


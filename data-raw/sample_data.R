library(vitality)

csv_fl_ls <- readRDS("../../grant_2026/data/clean/csv_fl_ls.rds")
taglife_rawDF <- csv_fl_ls$GPUD2026_taglife_17Aug2026

# Extract failure times by lot
lot1_data <- taglife_rawDF$days_difference[taglife_rawDF$lot == "Lot 1"]
lot2_data <- taglife_rawDF$days_difference[taglife_rawDF$lot == "Lot 2"]
lot3_data <- taglife_rawDF$days_difference[taglife_rawDF$lot == "Lot 3"]
pooled_data <- taglife_rawDF$days_difference


lot1_S <- fc_surv(time=lot1_data)#,rc.value =lot1_cens_val)
fc_plot(time = lot1_data, surv = lot1_S, hist=F, main="Lot 1")

fc_lot1_fit <- failCompare::fc_fit(time = sort(lot1_data), model = "vitality.ku", SEs = TRUE)
plot(fc_lot1_fit) # matches ATLAS
fc_lot1_fit


fc_lot1_fitRC <- failCompare::fc_fit(time = sort(lot1_data), model = "vitality.ku", SEs = TRUE,rc.value = 75)
plot(fc_lot1_fit) # matches ATLAS

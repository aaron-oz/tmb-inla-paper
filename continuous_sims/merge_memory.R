# this file merge memory use onto exp and iters
# it's a bit tricky b/c the task ids mask the exp and iter and because
# of the relaunches needed to achieve convergence
# aoz 2020

# dirs

mem.dir <- "~/Documents/Research/2020_tmb_v_inla/tmb_inla_sim/2020_08_08_13_13_30/"
map.dir <- "~/Documents/Research/2020_tmb_v_inla/tmb_inla_sim/2020_08_08_13_13_30/"

# set env
library(glue)
library(data.table)

# load memory
mem <- fread(file.path(mem.dir, "aaron-task-ram.csv"))

# load mapping
map <- fread(file.path(mem.dir, "jid_tid_TO_exp_iter_mapping.csv"))

# merge together
mem <- merge(mem, map, by.x = c("job_number", "task_number"), by.y = c("jid", "tid"), all.x = T, all.y = T)

# look at rows that have NAs, these all will be dropped
mem.na <- mem[which(apply(mem, 1, function(x){sum(is.na(x))}) >0),]
# not too bad - 17 rows dropped
mem <- na.omit(mem)

# ok, now we need to choose the LAST time a exp and iter show up. that is the final, aka successful run
# EXCEPT for those pesky NA rows. need to think about what to do with that...

# NOTE! no point in continuing. I realized that I only have the peak
# memory from the run. I don't kow if it was TMB or INLA that peaked
# it and I don't have any comparisons of INLA vs TMB...

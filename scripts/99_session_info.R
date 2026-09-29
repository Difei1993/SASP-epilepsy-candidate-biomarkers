# 99_session_info.R
# Record the R environment used for reproducibility.

source("scripts/00_config.R")

sink(file.path(RESULT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Wrote results/sessionInfo.txt")

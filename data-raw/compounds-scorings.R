compScorings <- read.csv(system.file("data-raw", "compounds-scorings.csv", package = "patRoon"), stringsAsFactors = FALSE)
usethis::use_data(compScorings, internal = TRUE, overwrite = TRUE)

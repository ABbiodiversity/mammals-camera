#-----------------------------------------------------------------------
# Packages
library(plumber)
library(here)

# Load model
model_moose <- readr::read_rds(here("R/moose_gap.rds"))

#* @apiTitle Probabilistic gap prediction

#* Give a prediction based on image gaps
#* @parser json
#* @post /gap_predict
function(req, res) {

  predict(model_moose, x = as.data.frame(req$body))

}

#* @plumber
function(pr) {

  pr %>%
    pr_set_api_spec(yaml::read_yaml(here("R/wildRtrax.yaml")))

}
















# JMP
Joint Model for Palliative care studies (R package)


## How to install the R package
require(devtools)

devtools::install_github("https://github.com/gitlzg/JMP", build_vignettes = TRUE)


## Detailed vignettes
browseVignettes("JMP")


## Example
library(JMP)

data(dat_JMP)

results <- JMP(dat=dat_JMP,
covariate.fields = c("sex", "qol_0"),
qol.prefix="qol_", qol.time.prefix="time_",
id.field = "id", survival.time.field = "survival_time",
censoring.status.field = "death", treatment.status.field = "trt")

## Extract monthly results from output
results$monthlyResults
## The overall p value for testing the entire model
results$overallPvalueModel


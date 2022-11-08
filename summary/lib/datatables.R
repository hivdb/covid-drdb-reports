library(DT)
library(slider)
library(dplyr)
library(kableExtra)

styledDT <- function(data, pageLength = 20) {
  data %>% datatable(
    escape = FALSE,
    rownames = FALSE,
    height = "auto",
    extensions = 'Buttons',
    options = list(pageLength = pageLength, dom = 'Bfrtip',
    buttons = list('copy', 'csv', list(extend = 'excel', title = NULL), 'print'))
  )
}
styledKBL <- function(data, caption = NA) {
  options(knitr.kable.NA = '')
  data %>% kbl(caption = caption, row.names = FALSE) %>% kable_minimal(full_width = FALSE)
}

simpleOrdinal <- function(x) {
  # since scales::ordinal doesn't handle decimals
  x = as.character(x)
  if (endsWith(x, "1")) {
    x = sprintf("%sst", x)
  }
  else if (endsWith(x, "2")) {
    x = sprintf("%snd", x)
  }
  else if (endsWith(x, "3")) {
    x = sprintf("%srd", x)
  }
  else {
    x = sprintf("%sth", x)
  }
  return(x)
}

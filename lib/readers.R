library(scales)

.getpwd <- function() {
  pwd = getSrcDirectory(.getpwd)[1]
  if (is.na(pwd) || startsWith(pwd, '.')) {
    pwd = '.'
  }
  print(pwd)
  return(pwd)
}

TABLES_DIR = normalizePath(sprintf('%s/../covid-drdb-payload/tables', .getpwd()))
LU_TABLES_DIR = normalizePath(sprintf('%s/../covid-drdb-payload/lu_tables', .getpwd()))

read.dbTable <- function(filename, colClasses = NA) {
  table_prefix = ifelse(
    startsWith(filename, "lu_"),
    yes = LU_TABLES_DIR,
    no = TABLES_DIR
  )
  read.table(
    sprintf("%s/%s", table_prefix, filename),
    header = TRUE,
    sep = ",",
    quote = '"',
    na.strings = "NULL",
    colClasses = colClasses
  )
}

read.dbTables <- function(dirname, colClasses = NA) {
  files = list.files(sprintf("%s/%s", TABLES_DIR, dirname),
                     pattern = "*.csv", full.names = FALSE)
  dfs = lapply(files, function(x) {
    read.dbTable(sprintf("%s/%s", dirname, x))
  })
  do.call(dplyr::bind_rows, dfs)
}

read.suscResults <- function(
  partialResistFold = 3,
  resistFold = 10
) {
  dfSusc = read.dbTables("susc_results", colClasses = c(fold = "numeric"))
  dfSusc$fold_cmp = ifelse(
    !is.na(dfSusc$fold) &
    !is.null(dfSusc$fold) & (
      dfSusc$fold_cmp == "" |
      is.na(dfSusc$fold_cmp)
    ),
    yes = "=",
    no = dfSusc$fold_cmp)
  dfSusc$ordinal_number = ifelse(
    is.na(dfSusc$ordinal_number),
    yes = 1,
    no = dfSusc$ordinal_number
  )
  dfSusc$cumulative_count = ifelse(
    is.na(dfSusc$cumulative_count),
    yes = 1,
    no = dfSusc$cumulative_count
  )
  dfSusc$resistance_level = ifelse(
    !is.na(dfSusc$resistance_level),
    dfSusc$resistance_level,
    ifelse(
      dfSusc$fold_cmp %in% c("<", "~", "=") &
      dfSusc$fold <= partialResistFold,
      "susceptible",
      ifelse(
        dfSusc$fold_cmp %in% c("<", "~", "=") &
        dfSusc$fold <= resistFold,
        "partial-resistance",
        "resistance"
      )
    )
  )
  dfSusc$resistance_level = factor(
    dfSusc$resistance_level,
    levels = c("susceptible", "partial-resistance", "resistance")
  )
  dfSusc
}

read.suscResultsCP <- function(
  partialResistFold = 3,
  resistFold = 10
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold
  )
  dfCP = read.dbTables("rx_conv_plasma")
  merge(dfSusc, dfCP, by = c("ref_name", "rx_name"))
}

read.suscResultsIP <- function(
  partialResistFold = 3,
  resistFold = 10
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold
  )
  dfIP = read.dbTables("rx_immu_plasma")
  merge(dfSusc, dfIP, by = c("ref_name", "rx_name"))
}

read.suscResultsMAb <- function(
  partialResistFold = 3,
  resistFold = 10
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold
  )
  dfMAb = read.dbTables("rx_antibodies")
  dfMAb = aggregate(ab_name ~ ref_name + rx_name, dfMAb, FUN = function(x) x)
  merge(dfSusc, dfMAb, by = c("ref_name", "rx_name"))
}
#-----------------------------------------------------------------------------------------------------------------------

# Title: Process GIS point summaries for each camera deployment
# Author: Marcus Becker

# Previous scripts:

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(DBI)
library(RSQLite)

# path to SC drive
sc <- "S:/GC_eric/FromEric/Sites_summaries/"

# Create connection
con_batch01 <- dbConnect(SQLite(), paste0(sc, "Round2020/20200224_SC_Sites_SummaryTables_Round_2020.sqlite"))
con_batch02 <- dbConnect(SQLite(), paste0(sc, "Round2020/20200429_SC_Sites_SummaryTables_Round_2020_batch02.sqlite"))
con_batch03 <- dbConnect(SQLite(), paste0(sc, "Round2020/20200506_SC_Sites_SummaryTables_Round_2020_batch03.sqlite"))

as.data.frame(dbListTables(con_batch02))

# Read table
pts <- dbReadTable(con_batch02, '20200429_OffGrid_Sites_SurveyYear_2019_2020_Points_batch02') # Takes a minute or two.



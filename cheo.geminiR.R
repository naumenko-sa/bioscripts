#using gemini database in R

setwd("~/Desktop/project_cheo/2016-10-28_gemini_test/")
library("RSQLite")

con = dbConnect(RSQLite::SQLite(),dbname="nu7823.db")
dbListTables(con)

variants = dbGetQuery(con,'select * from variants limit 1000')

dbClearResult(variants)
dbDisconnect(con)

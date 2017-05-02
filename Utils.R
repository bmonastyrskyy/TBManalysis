# Author : Bohdan Monastyrskyy
# Date : 2017-04-28
# Description : the file contains user-defined functions-utils 


# safe load packages
safe.load.package <- function(package.name){
  if (! require(package.name, character.only = TRUE)){
    install.packages(package.name);
    if (! require(package.name, character.only = TRUE)){
      stop(paste("Failed to load package", package.name));
    }
  }
}

# connect to database and return the connector in case of success
connect.database <- function(db.name, db.user, db.password, db.host, dp.port){
  require('RPostgreSQL')
  # loads the PostgreSQL driver
  drv <- dbDriver("PostgreSQL")
  # creates a connection to the postgres database
  con <- dbConnect(drv, dbname = db.name,
                   host = db.host, port = db.port,
                   user = db.user, password = db.password)
  # return
  con
}
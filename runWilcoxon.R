# Author : Bohdan Monastyrskyy
# Date : 2017-04-28
# Description : the script runs the Wilcoxon tests on composote z-scores
#         (gdt_ha , lddt, cad_aa, sg a.k.a sg_lvr_6_2, ace a.k.a qse)
#

# load config file, user-defined functions, and utils
source("config.R")
source("Utils.R")
source("Functions.R")

# load packages
LIBS <- c("RPostgreSQL", "dplyr")
tmp <- sapply(LIBS, safe.load.package)
rm(tmp)


# set up connection with database
con <- connect.database(db.name, db.user, db.password, db.host, dp.port)
rm(db.password) #delete db.password - security issue

# fetch data from database at predictioncenter.org
# df_1s <- fetchData(con, 'casp12', 'first', 'server')
df_1s <- fetchData(con, 'casp12', 'first', 'human')

# # add combined z-score
df_1s <- comp.z(df_1s)
#  add ranks
df_1s_ranks <- calcRanksGroups(df_1s)

m.w <- matrix('-', 20, 20)
m.t <- matrix('-', 20, 20)
m.h <- matrix('-', 20, 20)
m.dh <- matrix('-', 20, 20) # matrix of no_wins-no_losses
groups2skip <- c(48,342)
df_1s_ranks <- df_1s_ranks %>% filter(!gr_code %in% groups2skip)
top20 <- unlist(df_1s_ranks[1:20, c("gr_code")])
gr_names20 <- unlist(df_1s_ranks[1:20, c("gr_name")])

for (i in 1:min(20, length(top20))){
  cat(paste("\n"));
  for (j in 1:i){
    cat (paste(j," "))
    if (j < i) {
      wt <- my.wilcox(df_1s, top20[[i]], top20[[j]])
      if (is.null(wt)){
        m.w[j,i] <- sprintf('0')
      } else {
        m.w[j,i] <- sprintf("%d", wt$no.common)
        m.w[i,j] <- ifelse(wt$p.value < 0.01, sprintf("<0.01", wt$p.value), sprintf("%.3f", wt$p.value))
      }
      t <- my.t.test(df_1s, top20[[i]], top20[[j]])
      if (is.null(t)){
        m.t[j,i] <- sprintf('0')
      } else {
        m.t[j,i] <- sprintf("%d", t$no.common)
        m.t[i,j] <- ifelse(t$p.value < 0.01, sprintf("<0.01", t$p.value), sprintf("%.3f", t$p.value))
      }
      h2h <- head2head(df_1s, top20[[i]], top20[[j]])
      m.h[i,j] <- sprintf("%.3f", h2h$win)
      m.h[j,i] <- sprintf("%.3f", h2h$loss)
      m.dh[i,j] <- sprintf("%.3f", h2h$win - h2h$loss)
      m.dh[j,i] <- sprintf("%.3f", h2h$loss - h2h$win)
    }
  }
}
colnames(m.w) <- sapply(top20, FUN=function(x){sprintf("G%03d", x)})
rownames(m.w) <- mapply(FUN=function(x,y){sprintf("%s G%03d", x, y)}, gr_names20, top20 )

colnames(m.t) <- sapply(top20, FUN=function(x){sprintf("G%03d", x)})
rownames(m.t) <- mapply(FUN=function(x,y){sprintf("%s G%03d", x, y)}, gr_names20, top20 )


colnames(m.h) <- sapply(top20, FUN=function(x){sprintf("G%03d", x)})
rownames(m.h) <- mapply(FUN=function(x,y){sprintf("%s G%03d", x, y)}, gr_names20, top20 )

colnames(m.dh) <- sapply(top20, FUN=function(x){sprintf("G%03d", x)})
rownames(m.dh) <- mapply(FUN=function(x,y){sprintf("%s G%03d", x, y)}, gr_names20, top20 )

# write 
#write.csv(m.w, file = "server.top20.wilcox.csv", row.names = TRUE)
write.csv(m.w, file = "human.top20.wilcox.csv", row.names = TRUE)

#write.csv(m.t, file = "server.top20.ttest.csv", row.names = TRUE)
write.csv(m.t, file = "human.top20.ttest.csv", row.names = TRUE)

#write.csv(m.h, file = "server.top20.h2h.csv", row.names = TRUE)
write.csv(m.h, file = "human.top20.h2h.csv", row.names = TRUE)

#write.csv(m.dh, file = "server.top20.diff_h2h.csv", row.names = TRUE)
write.csv(m.dh, file = "human.top20.diff_h2h.csv", row.names = TRUE)


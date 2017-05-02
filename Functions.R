# Author : Bohdan Monastyrskyy
# Date : 2017-04-28
# Description : a collection of user-defined functions 
#     - to retrive and transform data
#     - carry out statistical tests
#     - generate summarized tables


fetchData <- function(con, casp='casp12',  bestfirst = 'first', humanserver='server' ){
  table <- 'zscores_1ms'
  if (humanserver == 'server'){
    if (bestfirst == 'first'){
      table <- 'zscores_1ms';
    } else {
      table <- 'zscores_bms';
    }
  } else {
    if (bestfirst == 'first') {
      table <- 'zscores_1m'
    } else {
      table <- 'zscores_bm'
    }
  }
  query <- paste0("select g.name as gr_name, g.code as gr_code, t.name || '-D' || z.domain as domain, dc.name as class, z.gdt_ha_4 as gdt_ha, z.lddt as lddt, z.cad_aa as cad, z.sg_lvr_6_2 as sg, z.qse as ase from ", casp, ".", table, " z ",
                  " join ", casp, ".groups g on g.id=z.groups_id ",
                  " join ", casp, ".domains d on z.targets_id=d.targets_id and d.index=z.domain ",
                  " join ", casp, ".domain_classifications dc on dc.id=d.domain_classifications_id ",
                  " join ", casp, ".targets t on t.id=d.targets_id ",
                  " where t.name similar to 'T0%' and t.cancellation_status=0 and d.index>0 and d.index<7 and dc.id in (1,2,4) ");
  df <- dbGetQuery(con, query);
  df
}

comp.z <- function(df){
  df$comp.z <- mapply(FUN=function(gdt_ha, lddt, cad, sg, ase) 
  {(0.333*gdt_ha + 0.111*lddt + 0.111*cad + 0.111*sg + 0.333*ase)},
  df$gdt_ha, df$lddt, df$cad, df$sg, df$ase)
  df
}

# rank groups according to formula
# z.gdt_ts +(z.lddt + z.cad + z.sg)/3 + z.ace
# all z-score < 0.0 are replaced by 0.0
calcRanksGroups <- function(df, classes=c("TBM", "FM/TBM")){
  tmp <- df
  tmp$comp.z0 <- mapply(FUN=function(gdt_ha, lddt, cad, sg, ase) 
  {(0.333*max(gdt_ha,0) + 0.111*max(lddt,0) + 0.111*max(cad,0) + 0.111*max(sg,0) + 0.333*max(ase,0))},
  tmp$gdt_ha, tmp$lddt, tmp$cad, tmp$sg, df$ase)
  tmp <- tmp %>% filter(class %in% classes) %>%  group_by(gr_name, gr_code) %>% summarize(avg_z0 = mean(comp.z0), no_dom=n(), sum_z0 = sum(comp.z0))
  tmp$rank <- rank(-tmp$sum_z0)
  tmp <- tmp[order(tmp$rank),]
  tmp
}

my.wilcox <- function(df, gr.code1, gr.code2, classes=c("TBM", "FM/TBM")){
  tmp1 <- df %>% filter(gr_code == gr.code1 & class %in% classes)
  tmp2 <- df %>% filter(gr_code == gr.code2 & class %in% classes)
  common.domains <- intersect(tmp1[, c("domain")], tmp2[, c("domain")])
  tmp1 <- tmp1 %>% filter( domain %in% common.domains)
  tmp1 <- tmp1[order(tmp1$domain),]
  tmp2 <- tmp2 %>% filter(domain %in% common.domains)
  tmp2 <- tmp2[order(tmp2$domain),]
  wt <- wilcox.test(unlist(tmp1$comp.z), unlist(tmp2$comp.z), paired = TRUE)
  wt$no.common <- length(common.domains)
  wt
}

my.t.test <- function(df, gr.code1, gr.code2, classes=c("TBM", "FM/TBM")){
  tmp1 <- df %>% filter(gr_code == gr.code1 & class %in% classes)
  tmp2 <- df %>% filter(gr_code == gr.code2 & class %in% classes)
  common.domains <- intersect(tmp1[, c("domain")], tmp2[, c("domain")])
  tmp1 <- tmp1 %>% filter( domain %in% common.domains)
  tmp1 <- tmp1[order(tmp1$domain),]
  tmp2 <- tmp2 %>% filter(domain %in% common.domains)
  tmp2 <- tmp2[order(tmp2$domain),]
  t <- t.test(tmp1$comp.z, tmp2$comp.z, paired = TRUE)
  t$no.common <- length(common.domains)
  t
}

head2head <- function(df, gr.code1, gr.code2, classes=c("TBM", "FM/TBM")){
  tmp1 <- df %>% filter(gr_code == gr.code1 & class %in% classes)
  tmp2 <- df %>% filter(gr_code == gr.code2 & class %in% classes)
  common.domains <- intersect(tmp1[, c("domain")], tmp2[, c("domain")])
  tmp1 <- tmp1 %>% filter( domain %in% common.domains)
  tmp1 <- tmp1[order(tmp1$domain),]
  tmp2 <- tmp2 %>% filter(domain %in% common.domains)
  tmp2 <- tmp2[order(tmp2$domain),]
  # number of wins of group gp.code1
  win1 <- 100.0*sum(unlist(mapply(FUN = function(x,y){x>y}, tmp1$comp.z, tmp2$comp.z)))/length(common.domains)
  # number of loss of group gp.code2
  loss1 <- 100.0*sum(unlist(mapply(FUN = function(x,y){x<y}, tmp1$comp.z, tmp2$comp.z)))/length(common.domains)
  list("win"= win1, "loss"=loss1)
}

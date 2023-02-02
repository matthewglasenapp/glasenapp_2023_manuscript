library(tidyverse)

# Code taken from Lanfear Lab http://www.robertlanfear.com/blog/files/concordance_factors.html

setwd(getwd())

d = read.delim("gcf.stat", header = T, comment.char='#')
#d = read.delim("scf.stat", header = T, comment.char='#')

chisq = function(DF1, DF2, N){
  tryCatch({
    # converts percentages to counts, runs chisq, gets pvalue
    chisq.test(c(round(DF1*N)/100, round(DF2*N)/100))$p.value
  },
  error = function(err) {
    # errors come if you give chisq two zeros
    # but here we're sure that there's no difference
    return(1.0)
  })
}               

e = d %>% 
  group_by(ID) %>%
  mutate(gEF_p = chisq(gDF1, gDF2, gN))
  #mutate(sEF_p = chisq(sDF1, sDF2, sN))

#subset(data.frame(e), (sEF_p < 0.05))

write.csv(e, "path", row.names=TRUE)

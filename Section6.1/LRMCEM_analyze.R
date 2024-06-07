library(dplyr);
library(reshape);

rm(list = ls())
df <- data.frame(read.table("summarytable.txt"));
names(df) <- c("iter", "eps", "epn", "b1", "b2", "b3", "tau", "m1", "m2");
head(df)

dfgroups <- df %>% group_by(epn) %>% summarise(size = n(), meanb1 = mean(b1), meanb2 = mean(b2), meanb3 = mean(b3), meantau = mean(tau),
                                               sdb1 = sd(b1), sdb2 = sd(b2), sdb3 = sd(b3), sdtau = sd(tau))
head(dfgroups)


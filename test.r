library(ltbt)

gamma.seq <- seq(-5, 5, 0.05)

exdata <- read.delim("RTEL1.tsv", sep="\t")

res <- run.ltbt(gamma.seq, exdata$gt, exdata$prs, exdata$outcome, 0.005)

res$h2Lx
# 0.1266585
res$pvalue
# 4.824949e-05
res$gamma.mle
# 1.45

plot(res$gamma, res$logL, ty="l", xlab="gamma", ylab="Log Likelihood Ratio")
abline(v=res$gamma.mle, col="red", lty=2)

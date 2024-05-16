


datamine <- data

test <- data[, c("ppid", "condition", "approach.avoid", "RT..s.")]
head(test)

colnames(test) <- c("subjects", "C", "R", "rt")


test$R <- factor(test$R)
test$subjects <- factor(test$subjects)
test$C <- factor(test$C)

head(test)
str(test)



dat <- test
fnams <- names(dat)[!(names(dat) %in% c("trials", "R", "rt"))]


a <- data.frame(table(test$R, test$C))
a$Freq <- a$Freq / sum(a$Freq) 

dens <- vector(mode = "list", length = dim(a)[1])
for (dims in 1:dim(a)[1]) {
  dens[[dims]] <- density(test[test$C == a[dims, 2] & test$R == a[dims, 1],]$rt, bw = "nrd0", adjust = 1)
  dens[[dims]]$y <- dens[[dims]]$y * a[dims, 3]
}
names(dens) <- c(paste(a$Var1, a$Var2, sep = "_"))


rx <- do.call(rbind, lapply(dens, function(x) {range(x$x)}) )
xlimi <- c(min(rx[, 1]), max(rx[, 2]))



densesy <- t(do.call(cbind, lapply(dens, function(x) {t(x[["y"]])})))
densesx <- t(do.call(cbind, lapply(dens, function(x) {t(x[["x"]])})))


resp <- rep(a$Var1, each = nrow(densesy) / nrow(a))
cond <- rep(a$Var2, each = nrow(densesy) / nrow(a))

d <- data.frame(densesy, densesx, resp, cond)
head(d)

library(ggplot2)

ggplot(d, aes(y =densesy, x = densesx, linetype = factor(resp), color = factor(cond))) +
  geom_line() + 
  theme_classic() +
  scale_linetype_manual(values = c("dashed", "solid"))

















a <- lapply(dens, function(x) {x[[y]]})
b <- a[[1]]

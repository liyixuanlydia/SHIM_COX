data(GBSG2,package="TH.data")
yv <- c(GBSG2$time,GBSG2$cens)
y <- matrix(yv,ncol=2,byrow=FALSE)
colnames(y) <- c("time","status")
GBSG2$tsizepnodes <- GBSG2$tsize * GBSG2$pnodes
GBSG2$tsizeestrec <- GBSG2$tsize * GBSG2$estrec
GBSG2$pnodesestrec <- GBSG2$pnodes * GBSG2$estrec
xdf <- data.frame(GBSG2$tsize,GBSG2$pnodes,GBSG2$estrec,GBSG2$tsizepnodes,GBSG2$tsizeestrec,GBSG2$pnodesestrec)
names(xdf) <- c("tsize","pnodes","estrec","tsize:pnodes","tsize:estrec","pnodes:estrec")
x <- as.matrix(xdf)
main_effect_names <- c("tsize","pnodes","estrec")
interaction_names <- c("tsize:pnodes","tsize:estrec","pnodes:estrec")

x <- x
y <- y
main.effect.names <- main_effect_names
interaction.names <- interaction_names
lambda <- c(0.0000005)
alpha <- c(0.5)
nlambda <- 1
threshold <- 1e-4
max.iter <- 2
center <- TRUE
normalize <- TRUE

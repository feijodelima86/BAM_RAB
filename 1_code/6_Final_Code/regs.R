#### P Plot ####

n1=2
n2=3
n3=1

par(mar=c(7, 7, 3, 3))

plot(P_TRD[,n1], P_TRD[,n2], 
     ylim = c(0, max(P_TRD[,n2],na.rm=TRUE)),
     xlim = c(min(P_TRD[,n1]), max(P_TRD[,n1],na.rm=TRUE)),
     cex=0,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     lwd=3
     
)

atx <- round(seq(min(P_TRD[,n1]), max(P_TRD[,n1], na.rm=TRUE), length.out=6), digits=2)
aty <- round(seq(0, max(P_TRD[,n2], na.rm=TRUE), length.out=6), digits=4)

Xlabel=expression(bold(paste("Turnover (%)")))
Ylabel=expression(bold("P Content (mg/g)"))

axis(side = 1, at = atx, labels=format(atx, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

title(xlab=Xlabel, line=2.5, cex.lab=1.5, family="Calibri")
title(ylab=Ylabel, line=2.7, cex.lab=1.5, family="Calibri")

pal <- colorRampPalette(c("#1b98e0", "red")) 

P_TRD$order = findInterval(P_TRD$P_TRD, sort(P_TRD$P_TRD))

points(P_TRD[,n1], P_TRD[,n2],pch=23, cex=2,TRD="black", bg=pal(nrow(P_TRD))[P_TRD$order],lwd=3)

box(lwd=2)

model <- lm(P_TRD[,n2] ~ P_TRD[,n1])

summary(model)

X <- range(P_TRD[,n1])
Y <- predict(model, newdata=data.frame(x=X))
lines(x=X, y=range(Y), lwd=2)

###########

#dev.new()

par(mar=c(7, 7, 3, 3))

plot(P_TRD[,n3], P_TRD[,n2], 
     ylim = c(0, max(P_TRD[,n2],na.rm=TRUE)),
     xlim = c(min(P_TRD[,n3]), max(P_TRD[,n3],na.rm=TRUE)),
     cex=0,
     xaxt = "n",
     xlab = NA,
     yaxt = "n",
     ylab = NA,
     lwd=3
     
)

atx <- round(seq(min(P_TRD[,n3]), max(P_TRD[,n3], na.rm=TRUE), length.out=6), digits=5)
aty <- round(seq(0, max(P_TRD[,n2], na.rm=TRUE), length.out=6), digits=4)

Xlabel=expression(bold(paste("P TRD (mg/L)")))
Ylabel=expression(bold("P Content (mg/g)"))

axis(side = 1, at = atx, labels=format(atx, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)
axis(side = 2, at = aty, labels=format(aty, scientific=F,digits = 1), las=1, font.axis=2, cex.axis=1)

title(xlab=Xlabel, line=2.5, cex.lab=1.5, family="Calibri")
title(ylab=Ylabel, line=2.7, cex.lab=1.5, family="Calibri")

pal <- colorRampPalette(c("#1b98e0", "red")) 

P_TRD$order = findInterval(P_TRD$TURNOVER, sort(P_TRD$TURNOVER))

points(P_TRD[,n3], P_TRD[,n2],pch=23, cex=2,TRD="black", bg=pal(nrow(P_TRD))[P_TRD$order],lwd=3)

box(lwd=2)

model <- lm(P_TRD[,n2] ~ P_TRD[,n3])

summary(model)

X <- range(P_TRD[,n3])
Y <- predict(model, newdata=data.frame(x=X))
lines(x=X, y=range(Y), lwd=2)



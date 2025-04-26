
pdf("Beta_Dens_Symm.pdf", height = 8, width = 12)
curve(dbeta(x, 1, 1), col=1, ylim=c(0,3), ylab="Density", lwd=1.5)
curve(dbeta(x, 2, 2), col=2, add=T, lwd=1.5)
curve(dbeta(x, 5, 5), col=4, add=T, lwd=1.5)
curve(dbeta(x, .5, .5), col=3, add=T, lwd=1.5)
legend(.6,3, c("a = b = 1",
               "a = b = 2",
               "a = b = 5",
               "a = b = 0.5"),
       col=c(1,2,4,3), lty=1, lwd=2, bty = "n", cex = 1.5)
dev.off()

pdf("Beta_Dens_Asymm.pdf", height = 8, width = 12)
curve(dbeta(x, 1, 1), col=1, ylim=c(0,3), ylab="Density", lwd=1.5)
curve(dbeta(x, 5, 2), col=2, add=T, lwd=1.5)
curve(dbeta(x, .5, 2), col=3, add=T, lwd=1.5)
legend(.6,3, c("a = b = 1",
               "a = 5, b = 2",
               "a = 0.5 = 2"),
       col=c(1,2,3), lty=1, lwd=2, bty = "n", cex = 1.5)
dev.off()

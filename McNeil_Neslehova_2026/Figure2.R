library(udp)

pdf("TplotsLegendre.pdf",width=8, height=2)
par(mfrow = c(1,5),     #
    oma = c(1, 1, 1, 1), # two rows of text at the outer left and bottom margin
    mar = c(2, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(1.2, 0.5, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA, pty ="s")
for (j in 2:6)
  plot(udplegendre(j), embellish = "bw")
par(mfrow = c(1,1))
dev.off()

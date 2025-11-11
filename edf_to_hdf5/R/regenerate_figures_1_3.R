# file regenerate_figures_1_3.R
# copyright (C) 2022-2025 Artur Gramacki and Jarosław Gramacki
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

source("funs.R")
source("my_custom_palette.R")


#///////////////////////////////////////////////////////////////////////////////////////////
# Figure 1 ----
#///////////////////////////////////////////////////////////////////////////////////////////
N <- 512
fs <- 256

par(mfrow = c(2,2), pty = "m")
par(mai = c(0.4, 0.4, 0.3, 0.2))

n <- 4

sigmas <- c(0.5, 0.2, 0.8, 0.5)
frequencies <- c(14, 8, 4, 1)
phases <- c(0, 1, 1.5, -2)
means = c(0.5, 0.8, 1, 1.5)

for (i in c(1:n)) {
  sigma <- sigmas[i]
  frequency <- frequencies[i]
  phase <- phases[i]
  mean <- means[i]
  main <- TeX(paste(
    "$\\mu=$", means[i], ", ", 
    "$\\sigma=$", sigmas[i], ", ",
    "$\\f=$", frequencies[i], ", ",
    "$\\phi=$", phases[i], 
    sep = ""
  )
  )
  
  gb <- Gabor.fun(N, fs, mean, phase, sigma, frequency, normalization = F)
  crossprod(gb$gabor)
  plot(
    gb$t, 
    gb$gauss, type = "l", ylim = c(-1, 1), col = "grey", 
    xlab = "", ylab = "",
    xaxt = "t", yaxt = "t", bty = "or",
    cex.axis = 1,
    main = main)
  lines(gb$t, gb$cosinus, type = "l", col = "grey")
  lines(gb$t, gb$gabor, col="black", lwd = 2)
}


#///////////////////////////////////////////////////////////////////////////////////////////
# Figure 2a ----
#///////////////////////////////////////////////////////////////////////////////////////////
SQLiteFile = "sample/data_for_figures_2.db"

out <- mp2tf(
  SQLiteFile, 
  channel = 12, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 8,
  displayCrosses = TRUE, 
  grid = FALSE, 
  crossesColor = "white", 
  palette = "my custom palette",
  rev = T,
  outMode = "plot",
  plotSignals = FALSE
)

#///////////////////////////////////////////////////////////////////////////////////////////
# Figure 2b ----
#///////////////////////////////////////////////////////////////////////////////////////////
out2 <- mp2tf(
  SQLiteFile, 
  channel = 12, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 16,
  outMode = "RData",  
  fileName = "sample/Figure_2.RData",
  fileSize = c(64, 64)
)

load("sample/Figure_2.RData")
par(mfcol = c(1, 1), pty = "m")
par(pty = "s")
par(mai = c(0.9, 0.9, 0.2, 0.4))
image(
  x = seq(0, 10, length.out = 64), 
  y = seq(0, 32, length.out = 64), 
  z = tf.matrix, 
  col = my_custom_palette,
  xlab = "Time [s]",
  ylab = "Frequency [Hz]",
)

file.remove("sample/Figure_2.RData")

#///////////////////////////////////////////////////////////////////////////////////////////
# Figure 3 ----
#///////////////////////////////////////////////////////////////////////////////////////////
par(mfcol = c(8, 1), pty = "m")
par(mai = c(0, 0.8, 0, 0))
par(mgp = c(0, 1, 0))

range <- range(out$originalSignal)

atoms <- matrix(nrow = nrow(out$gaborFunctions), ncol = out$epochSize)
for (i in 1:50) {
  atoms[i,] <- out$gaborFunctions[i,] * sqrt(out$atoms$energy[i] * out$f)  
}

plot(out$originalSignal, type = "l", xaxt = "n", yaxt = "n", bty = "n",  xlab = "n", ylab = "")
mtext("Signal", side = 2, line = 0, las = 2)
plot(out$reconstructedSignal, type = "l", xaxt = "n", yaxt = "n", bty = "n",  xlab = "", ylab = "")
mtext("Reconstr.", side = 2, line = 0, las = 2)
plot(atoms[1,], type = "l", ylim = range, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", col = "blue")
mtext("Atom 1", side = 2, line = 0, las = 2)
plot(atoms[2,], type = "l", ylim = range, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", col = "blue")
mtext("Atom 2", side = 2, line = 0, las = 2)
plot(atoms[3,], type = "l", ylim = range, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", col = "blue")
mtext("Atom 3", side = 2, line = 0, las = 2)
plot(atoms[4,], type = "l", ylim = range, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", col = "blue")
mtext("Atom 4", side = 2, line = 0, las = 2)
plot(atoms[5,], type = "l", ylim = range, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", col = "blue")
mtext("Atom 5", side = 2, line = 0, las = 2)

sum <- 0
for (j in 1:20) {
  sum <- sum + atoms[j,]
}
plot(sum, type = "l", ylim = range, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "" , col = "red")
mtext("A1-A20", side = 2, line = 0, las = 2)


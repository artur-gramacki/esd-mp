# file bin_to_RData_STFT.R
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

patterns <- c(".*fs_eA.*.bin", ".*fns_eA.*.bin", ".*fs_eB.*.bin", ".*fns_eB.*.bin", ".*fs_eC.*.bin", ".*fns_eC.*.bin")
dirs <- c("eA_w10_c20", "eA_w10_c20", "eB_w10_c20", "eB_w10_c20", "eC_w10_c20", "eC_w10_c20")  

dir <- "../"

save_to_png <- FALSE

if (save_to_png) {
  ext <- ".png"
  } else { ext <- "png"
  ext <- ".RData"
  }

tf_file_size = c(64, 64)

fmax <- 32

for (i in 1:6) {
  cat(dirs[i], "\n", sep = "")
  files <- list.files(paste(dir, "working/bin_files/", dirs[i], sep = ""), include.dirs = FALSE, pattern <- patterns[i])
  # remove the extension
  files <- tools::file_path_sans_ext(files)

  for (j in 1:length(files)) {
    cat(j, " / ", length(files), "\n", sep = "")
    fn <- paste(dir, "working/bin_files/", dirs[i], "/", files[j], ".bin", sep = "")
    fs <- file.size(fn)
    data <- NA

    fh = file(fn, "rb")
    for (k in 1:(fs / 8)) {# 8: size of double in bytes
      (b <- readBin(con = fh, what = "raw", n = 8)) 
      # swap to use little-endian
      (hex_str <- paste(b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], sep = ""))
      # https://stackoverflow.com/questions/39461349/converting-hex-format-to-float-numbers-in-r
      
      # Convert to raw
      (raw_vec <- as.raw(strtoi(substring(hex_str, (step <- seq(1, nchar(hex_str), by = 2)), step + 1), 16)))
      # Convert raw to double
      (num <- readBin(raw_vec, what = "double", n = 1, size = 8, , endian = "little"))
      data[k] <- num
    }
    
    close(fh)
    data.mtx <- matrix(data, nrow = fs / 8 / 18, ncol = 18, byrow = TRUE)

    for (m in 1: 18) { # we have 18 channels in our files
      tfFileName <- paste(
        dir,
        "working/tf_maps_STFT/",
        dirs[i],
        "_",
        tf_file_size[1],
        "_",
        tf_file_size[2],
        "/",
        files[j],
        "_ch",
        formatC(m, width = 2, format = "d", flag = "0"),
        "_",
        tf_file_size[1],
        "_",
        tf_file_size[2],
        "_STFT",
        "_file",
        formatC(j, width = 5, format = "d", flag = "0"),
        ext, 
        sep = "")
      # -- for testing only --
      # cat (tfFileName, "\n")
      # plot(data.mtx[,m], type ="l", main = tfFileName)
      
      # It seems that e1071::stft works better than gsignal::stft
      # ft <- gsignal::stft(x = data.mtx[,m], fs = 256, overlap = 0.99, window = 128)

      # Much slower but gives t-f maps with better resolution in f axis
      #ft <- e1071::stft(X = data.mtx[,m], win = 128, inc = 1, coef = 512)
      ft <- e1071::stft(X = data.mtx[,m], win = 128, inc = 1, coef = 64)
      
      # y axis limit to fmax Hz, our EEG recordings' sampling freq = 256
      # 256 / 2: Nyquist freq
      p <- (ncol(ft$values) * fmax) / (256 / 2)

      if (save_to_png) {
        graphics.off()
        png(tfFileName, width = 640, height = 640, pointsize = 18)
        #par(pty = "m", mai = c(0, 0, 0, 0))

        graphics::image(
          z = ft$values[,1:p], 
          col = my_custom_palette,
          xlab = "Time [s]",
          ylab = "Frequency [Hz]",
          xaxt = "n",
          yaxt = "n",
          main = tfFileName,
          cex.main = 0.8
        )
        axis(
          side = 1, 
          at = seq(from = 0, to = 1, length.out = 5), 
          labels = seq(from = 0, to = 10, length.out = 5)
        )
        axis(
          side = 2, 
          at = seq(from = 0, to = 1, length.out = 9), 
          labels = seq(from = 0, to = 32, length.out = 9)
        )
        dev.off()
      } # if (save_to_png)
      
      if (!save_to_png) {
        zz <- ft$values[,1:p]      
        rr <- raster::raster(nrow = ncol(zz), ncol = nrow(zz)) # this is how it should be: nrow = ncol(zz), ncol = nrow(zz)
        rr[] <- t(zz)
        tt <- raster::raster(ncol = tf_file_size[1], nrow = tf_file_size[2])
        tt <- raster::resample(rr, tt)
        m2 <- matrix(tt@data@values, tf_file_size[1], tf_file_size[2])
        # -- for testing only --
        # graphics::image(m2, col = my_custom_palette)
        
        # Rescaling to the range 0-1
        # Protect against a situation where a zero appears in the denominator
        if (max(m2) - min(m2) == 0) {
          tf.matrix <- matrix(0, tf_file_size[1], tf_file_size[2])  
        } else {
          tf.matrix <- (m2 - min(m2)) / (max(m2) - min(m2))
        }
        save(tf.matrix, file = tfFileName)
      
        # -- for testing only --
        # load(tfFileName)
        # graphics::image(tf.matrix, col = my_custom_palette)
      } # if (!save_to_png)
    } # for (m in 1: 18)

    df <- data.frame(
       i = paste(i, " / 6", sep = ""),
       j = paste(j, " / ", length(files), sep = ""),
       row.names = NULL)
    write.table(
      df, 
      paste(dir, "working/tf_maps_STFT/progress.txt", sep = ""),
      row.names = F, 
      col.name = TRUE, 
      sep = "\t", 
      quote = FALSE)
    
  } # for (j in 1:length(files))
} # for (i in 1:3) 




# file funs.R
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

# These two libraries must by first installed.
# Install from https://bioconductor.org
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("rhdf5")
to_install = !pkgs %in% installed.packages()
if(to_install) {
  BiocManager::install("rhdf5")
}

pkgs <- c(
  "grid", "hht", "RSQLite", "grDevices", "colorRamps", 
  "gplots", "fields", "raster", "DescTools", "mvtnorm",
  "signal", "Matrix", "latex2exp", "edf", "png")

to_install = !pkgs %in% installed.packages()
if(any(to_install)) {
  install.packages(pkgs[to_install])
}

library('grid')
library('hht') 
library("RSQLite")
library("grDevices") 
library("colorRamps") 
library("gplots") 
library("fields") 
library("raster") 
library("DescTools") 
library("mvtnorm")
library("signal") 
library('Matrix') 
library("latex2exp")
library("edf")
library("png")
library("rhdf5")
library("signal")

# ///////////////////////////////////////////////////////////////////////////////////////////
select_seizure_chunks <- function(data, f, k) {
  # data - seizure annotations file
  # f - base frequency 
  # k - which patient (1:79 for our EEG database)
  
  aa <- 0
  seizures <- data.frame()
  sec.1 <- which(data[,k] == 1)
  # only for patients with at least one seizure 
  if (length(sec.1) != 0 ) {
    # https://stackoverflow.com/questions/23095415/how-to-find-if-the-numbers-are-continuous-in-r
    # s: each list item is the starting and ending second of an epileptic seizure
    s <- unname(tapply(sec.1, cumsum(c(1, diff(sec.1)) != 1), range))
    for (i in 1:length(s)) {
      s.duration.in.secs <- s[[i]][2] - s[[i]][1] + 1
      s.samples <- c( 
        (s[[i]][1] - 1) * f + 1, 
        (s[[i]][1] - 1)* f  + (s.duration.in.secs * f)
      )
      seizures[i, 1] <- k
      seizures[i, 2] <- s.duration.in.secs
      seizures[i, 3] <- s[[i]][1]
      seizures[i, 4] <- s[[i]][2]
      seizures[i, 5] <- s.samples[1]
      seizures[i, 6] <- s.samples[2]
    }  
    colnames(seizures) <- c("patient", "seizure_duration", "from_sec", "to_sec", "from_sample", "to_sample")
    list(
      seizures = seizures, 
      total.seizures = length(s)
    )
  } else {
    list(seizures = NA, total.seizures = NA)
  }  
}  

# ///////////////////////////////////////////////////////////////////////////////////////////
generate_montage <- function(matrix) {
  # The code based on read_data_montage.m  
  # https://github.com/ktapani/Neonatal_Seizure_Detection
  str <- data.frame()
  str[1,1] = 'Fp2'; str[1,2] = 'F4';    # Fp2-F4
  str[2,1] = 'F4'; str[2,2] = 'C4';     # F4-C4
  str[3,1] = 'C4'; str[3,2] = 'P4';     # C4-P4
  str[4,1] = 'P4'; str[4,2] = 'O2';     # P4-O2
  str[5,1] = 'Fp1'; str[5,2] = 'F3';    # Fp1-F3
  str[6,1] = 'F3'; str[6,2] = 'C3';     # F3-C3
  str[7,1] = 'C3'; str[7,2] = 'P3';     # C3-P3
  str[8,1] = 'P3'; str[8,2] = 'O1';     # P3-O1
  str[9,1] = 'Fp2'; str[9,2] = 'F8';    # Fp2-F8
  str[10,1] = 'F8'; str[10,2] = 'T4';   # F8-T4
  str[11,1] = 'T4'; str[11,2] = 'T6';   # T4-T6
  str[12,1] = 'T6'; str[12,2] = 'O2';   # T6-O2
  str[13,1] = 'Fp1'; str[13,2] = 'F7';  # Fp1-F7
  str[14,1] = 'F7'; str[14,2] = 'T3';   # F7-T3
  str[15,1] = 'T3'; str[15,2] = 'T5';   # T3-T5
  str[16,1] = 'T5'; str[16,2] = 'O1';   # T5-O1
  str[17,1] = 'Fz'; str[17,2] = 'Cz';   # Fz-Cz
  str[18,1] = 'Cz';  str[18,2] ='Pz';   # Cz-Pz
  
  label <- colnames(matrix)
  mtx.mont <- matrix(NA, nrow = nrow(matrix), ncol = 20)
  
  for (jj in 1:18) {
    ref1 = rep(0, 21)
    ref2 = rep(0, 21)
    for (ii in 1:21) {
      ref1[ii] <- as.numeric(grepl(str[jj,1], label[ii]))
      ref2[ii] <- as.numeric(grepl(str[jj,2], label[ii]))
    }
    qq1 = which(ref1 == 1)[1]
    qq2 = which(ref2 == 1)[1]
    mtx.mont[,jj] = matrix[, qq1] - matrix[, qq2] 
  }
  df.mont <- as.data.frame(mtx.mont)
  df.mont[,19] = matrix[,22]
  for (kk in 1:18) {
    colnames(df.mont)[kk] <- paste(str[kk,1], "-", str[kk,2], sep ="")
  }
  colnames(df.mont)[19] <- c("t")
  df.mont
}

# ///////////////////////////////////////////////////////////////////////////////////////////
filters_coeff <- function(fs = 256, notch = c(48.5, 51.5), lowpass = 30, highpass = 1) {
  # https://openbci.com/forum/index.php?p=/discussion/193/50hz-notch-filter-coefficients

  ## 50 Hz notch filter
  bf.notch <- butter(2, notch / (fs / 2), "stop")
  freqz(bf.notch)

  # Low pass IIR Butterworth, cutoff at 'lowpass' Hz
  bf.low <- butter(8, lowpass / (fs / 2), "low")
  freqz(bf.low)

  
  # High pass IIR Butterwoth, cutoff at 'highpass' Hz
  bf.high <- butter(2, highpass / (fs / 2), "high")
  freqz(bf.high)

  
  list(bf.notch = bf.notch, bf.low = bf.low, bf.high = bf.high)
}

# ///////////////////////////////////////////////////////////////////////////////////////////
generate_samples <- function(
  which.expert,
  annotations_file,
  seizure.IDs,
  non.seizure.IDs,
  window,
  chunks,
  down.sampling.factor = 1,
  filtering = FALSE,
  dir,
  random = FALSE,
  my.seed = 42,
  write.txt.files = FALSE,
  write.hdf5.files = FALSE,
  write.bin.files = FALSE,
  write.bin.as.txt.files = FALSE) {
  
  # ///////////////////////////////////////////////////////////////////////////////////////////
  # which.expert - "A", "B" or "C" (or any other symbol) 
  #
  # annotations_file - file name where annotations are saved
  #
  # seizure.IDs - infant IDs which have seizures. 
  #               We assume that EDF file names are eeg1.edf, eeg2.edf, eeg3.edf etc.
  #
  # non.seizure.IDs - infant IDs which are seizure free. 
  #                   We assume that EDF file names are eeg1.edf, eeg2.edf, eeg3.edf etc.
  #
  # window - window size in seconds
  #
  # chunks - number of window's chunks  
  #
  # down.sampling.factor - down sampling factor of the original edf file   
  #                        (must be divisible by "f", that is, 1,2,4,8 etc. 
  #                        In practice no more than 4, 1 - no down-sampling)
  #
  # filtering - see preprocess.m in https://github.com/ktapani/Neonatal_Seizure_Detection
  #
  # dir - root dir, see the directory structure given in the paper
  #
  # random - if TRUE generate random non-seizure chunks from non-seizure EDF files 
  #
  # my.seed - for selecting non-seizures chunks
  #
  # write.txt.files - if TRUE safe the seizure and nos-seizure final files as TXT files
  #                   (we generate TXT files for illustrative purposes only. As for the content, 
  #                   they are fully compatible with HDF5 binary files)     
  #
  # write.hdf5.files - if TRUE safe the file in the HDF5 format
  #
  # write.bin.files - if TRUE save binary files (doubles) in the format required by empi
  # 
  # write.bin.as.text.files - if TRUE save binary files also as text files 
  #                           (only for controlling purposes)
  #
  # /////////////////////////////////////////////////////////////////////////////////////////// 
  
   fout <- filters_coeff(fs = 256, notch = c(48.5, 51.5), lowpass = 30, highpass = 1)
   bf.notch = fout$bf.notch
   bf.low = fout$bf.low
   bf.high = fout$bf.high 
   
  # read annotations file  
  ann <-
    read.csv(
      paste(dir, "annotations/", annotations_file, sep = ""),
      sep = ",",
      header = TRUE,
      stringsAsFactors = F,
      check.names = FALSE,
      encoding = 'UTF-8'
    )
  
  # Read EDF params from the first file.
  # n.sigs - number of recorded channels in the EDF file (+1 because then we add "t" column) 
  # f.edf - sampling frequency of the original EDF files
  files <- list.files(paste(dir, "edf", sep = ""))
  filename <- paste(dir, "edf/", files[1], sep = "")
  edf <- read.edf(filename = filename, read.annotations = FALSE, header.only = FALSE)
  n.sigs <- edf[["header.global"]][["n.signals"]] + 1
  f.edf <- edf[["header.signal"]][[1]][["n.samples"]]
  f <- f.edf /down.sampling.factor
  
  cat("---------------------------------------------------------------------------", "\n")
  cat("annotations file name:     ", annotations_file, "\n", sep = "")
  cat("seizure patients:         ", seizure.IDs, "\n", sep = " ")
  cat("non-seizure patients:     ", non.seizure.IDs, "\n", sep = " ")    
  if (write.hdf5.files) cat("hdf5 file being generated: ", "expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.hdf5", "\n", sep = "")
  if (write.txt.files) cat("txt file being generated:  ", "expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", "\n", sep = "")
  cat("---------------------------------------------------------------------------", "\n\n")
  
  SEIZURE <- data.frame()
  NON.SEIZURE <- data.frame() 
  FINAL <- data.frame()
  
  unlink(paste(dir, "working/aux_files/seizures_",     which.expert, "_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/aux_files/seizures_",     which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/aux_files/non_seizures_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  
  unlink(paste(dir, "working/aux_files/SEIZURE_expert_",      which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/aux_files/NON.SEIZURE_expert_",  which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))
  unlink(paste(dir, "working/aux_files/expert_",              which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""))

  # Uncomment if you want to delete all files before generating new ones
  # files <- list.files(paste(dir, "working/bin_files/", sep = ""), include.dirs = TRUE)
  # if (length(files) != 0) file.remove(paste(dir, "working/bin_files/", files, sep = ""))
  
  # Counter of seizure chunks. 
  S <- 0
  
  # ///////////////////////////////////////////////////////////////////////////////////////////
  # First, we analyze patients with at least one annotated seizure
  # ///////////////////////////////////////////////////////////////////////////////////////////
  for (i in seizure.IDs) {
    filename <- paste(dir, "edf/eeg", i, ".edf", sep = "")
    cat(filename, sep = "")
    edf <- read.edf(filename = filename, read.annotations = FALSE, header.only = FALSE)
    
    # Calculate the number of samples. Take the first signal because each has the same number of samples. 
    len <- length(edf$signal[[1]]$data)
    
    edf.mtx <- matrix(NaN, nrow = len, ncol = n.sigs)
    for (n in 1:(n.sigs - 1)) {
      edf.mtx[, n] <- edf[["signal"]][[n]][["data"]]
    }
    head(edf.mtx)
    
    sig.names <- NaN
    # Read the names of the signals.
    for (s in 1:(n.sigs - 1)) {
      sig.names[s] <- edf$header.signal[[s]]$label
    }
    colnames(edf.mtx) <- c(sig.names, "t")
    head(edf.mtx)
    
    # Add time stamps to the last column. Take data from any channel, the same everywhere. 
    edf.mtx[, n.sigs] <- edf[["signal"]][[1]][["t"]]
    head(edf.mtx)
    
    # down-sampling
    edf.mtx <- edf.mtx[seq(1, nrow(edf.mtx), down.sampling.factor),]    
    
    # montage
    edf.mtx.m <- generate_montage(edf.mtx)
    # last column - reserved for class label
    # penultimate column - timestamps
    head(edf.mtx.m)
    cc <- ncol(edf.mtx.m)
    
    # ing, see preprocess.m in https://github.com/ktapani/Neonatal_Seizure_Detection  
    if (filtering) {
       edf.mtx.p <- matrix(NA, nrow = nrow(edf.mtx.m), ncol = ncol(edf.mtx.m))
      colnames(edf.mtx.p) <- colnames(edf.mtx.m)
      for (m in 1:(cc - 2)) {
        edf.mtx.p[, m] = signal::filtfilt(bf.notch, edf.mtx.m[, m]); # 50Hz notch filter
        edf.mtx.p[, m] = signal::filtfilt(bf.low, edf.mtx.p[, m]); # Low pass IIR Butterworth
        edf.mtx.p[, m] = signal::filtfilt(bf.high, edf.mtx.p[, m]); # High pass IIR Butterwoth
      }
      edf.mtx.p[, cc - 1] <- edf.mtx.m[, cc - 1]
    } else {
      edf.mtx.p <- edf.mtx.m
    }
    head(edf.mtx.p)
    
    out <- select_seizure_chunks(data = ann, f = f, k = i)
    out$seizures
    out$total.seizures
    cat(" (total number of seizures annotated: ", out$total.seizures, ")\n", sep = "")
    
    if (i == seizure.IDs[1]) cn = TRUE else cn = FALSE
    
    # We are temporarily disabling warnings to avoid displaying an unobtrusive warning
    options(warn = -1)
    
    write.table(
      x = out$seizures,
      file = paste(dir, "working/aux_files/seizures_", which.expert, "_", f, "Hz.txt", sep = ""),
      append = TRUE,
      col.names = cn,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    options(warn = 0)
    
    seizure.mtx <- matrix(NaN, nrow = 0, ncol = ncol(edf.mtx.p)) 
    
    m <- 0
    out2 <- data.frame()
    for (j in 1:out$total.seizures) {
      (r <- out$seizures[j, 2] %/% window)
      if (r > 0) {
        if (r > chunks) { # so as not to select more chunks than possible 
          fr <- chunks
        } else {
          fr <- r
        }
        for (k in 1:fr) {
          S <- S + 1
          m <- m + 1
          (from <- out$seizures[j, 5] + (k - 1) * f * window) 
          (to <- out$seizures[j ,5] + (k * window * f) - 1)  
          out2[m, 1] <- out$seizures[j, 1]
          out2[m, 2] <- out$seizures[j, 2]
          out2[m, 3] <- from
          out2[m, 4] <- to
          out2[m, 5] <- k
          out2[m, 6] <- m
        }  
      } else {
        # do nothing
      }
    }
    colnames(out2) <- c("patient", "seizure_duration", "from_sample",	"to_sample",  "chunk#", "seq")
    
    if (i == seizure.IDs[1]) cn = TRUE else cn = FALSE
    
    # We are temporarily disabling warnings to avoid displaying an unobtrusive warning
    options(warn = -1)
    
    write.table(
      x = out2,
      file = paste(dir, "working/aux_files/seizures_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      append = TRUE,
      col.names = cn,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    options(warn = 0)
    
    if (nrow(out2) > 0) {
      for (p in 1:nrow(out2)) {
        (from <- out2[p, 3])
        (to <- out2[p, 4])
        (chunk <- out2[p, 5])
        (seq <- out2[p, 6])
        temp <- edf.mtx.p[from:to,]
        
        if (filtering) txt = "f" else txt = ""
        ext <- ".txt"
        file = paste(
          dir, 
          "working/bin_files/",
          txt,
          "s", 
          "_e", which.expert, 
          "_p", formatC(i, width = 2, format = "d", flag = "0"), 
          "_w",  formatC(window, width = 2, format = "d", flag = "0"),  
          "c", formatC(chunks, width = 2, format = "d", flag = "0"),
          "_c", formatC(chunk, width = 2, format = "d", flag = "0"), 
          "_seq_", formatC(seq, width = 4, format = "d", flag = "0"),
          "_",
          f, "Hz", sep = "") 
        file.txt <- paste(file, ext, sep = "")
        
        # Write to a text file (just to make it easier to view the content, we don't really need this file)
        if (write.bin.as.txt.files) {
          write.table(
            x = round(temp, 4),
            file = file.txt,
            dec = ".",
             sep = "\t",
             col.names = TRUE,
             row.names = FALSE
          )
        }

        # To binary file
        # https://gregstoll.com/~gregstoll/floattohex/
        if (write.bin.files) {
          ext = ".bin"
          file.bin <- paste(file, ext, sep = "")
          fh = file(file.bin, "wb")
          for (n in 1:nrow(temp)) {
            # "-2", because the last two columns contain the timestamp and class label
            data_row <- as.numeric(temp[n, 1:(ncol(edf.mtx.p) - 2)])
            writeBin(data_row, fh) 
          }
          close(fh)
        }
        
        seizure.mtx <- rbind(seizure.mtx, temp)
      }
      head(seizure.mtx)
      nrow(seizure.mtx)
      
      seizure.mtx[, cc] <- rep(1, nrow(seizure.mtx))
      colnames(seizure.mtx)[cc] <- "seizure"
      
      SEIZURE <- rbind(SEIZURE, seizure.mtx)
    }
    
  } ### for (i in seizure.IDs)
  
  if (write.txt.files) {
    write.table(
      x = round(SEIZURE, 4),
      file = paste(dir, "working/aux_files/SEIZURE_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      dec = ".",
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
  }
  
  # ///////////////////////////////////////////////////////////////////////////////////////////
  # Second, we analyze patients with NO ONE annotated seizure. 
  # We select RANDOMLY as many chunks as we have selected from the seizured files. 
  # ///////////////////////////////////////////////////////////////////////////////////////////
  # Set how many chunks to take from each non-seizured file, so that there are about the same number of chunks.
  nn <- ceiling(S / length(non.seizure.IDs))

  for (q in non.seizure.IDs) {
    (filename <- paste(dir, "edf/eeg", q, ".edf", sep = ""))
    cat(filename, "\n", sep = "")
    edf <- read.edf(filename = filename, read.annotations = FALSE, header.only = FALSE)
    
    # Calculate the number of samples. Take the first signal because each has the same number of samples. 
    len <- length(edf$signal[[1]]$data)

    edf.mtx <- matrix(NaN, nrow = len, ncol = n.sigs)
    for (n in 1:(n.sigs - 1)) {
      edf.mtx[, n] <- edf[["signal"]][[n]][["data"]]
    }
    head(edf.mtx)

    sig.names <- NaN
    # Read the names of the signals.
    for (s in 1:(n.sigs - 1)) {
      sig.names[s] <- edf$header.signal[[s]]$label
    }
    colnames(edf.mtx) <- c(sig.names, "t")
    head(edf.mtx)

    # Add time stamps to the last column. Take data from any channel, the same everywhere. 
    edf.mtx[, n.sigs] <- edf[["signal"]][[1]][["t"]]
    head(edf.mtx)

    # down-sampling 
    edf.mtx <- edf.mtx[seq(1, nrow(edf.mtx), down.sampling.factor),]

    # montage
    edf.mtx.m <- generate_montage(edf.mtx)
    head(edf.mtx.m)
    cc <- ncol(edf.mtx.m)

    # filtering  
    if (filtering) {
      edf.mtx.p <- matrix(NA, nrow = nrow(edf.mtx.m), ncol = ncol(edf.mtx.m))
      colnames(edf.mtx.p) <- colnames(edf.mtx.m)
      for (m in 1:(cc - 2)) {
        edf.mtx.p[, m] = signal::filtfilt(bf.notch, edf.mtx.m[, m]); # 50Hz notch filter
        edf.mtx.p[, m] = signal::filtfilt(bf.low, edf.mtx.p[, m]); # Low pass IIR Butterworth
        edf.mtx.p[, m] = signal::filtfilt(bf.high, edf.mtx.p[, m]); # High pass IIR Butterwoth
      }
      edf.mtx.p[, cc - 1] <- edf.mtx.m[, cc - 1]
    } else {
      edf.mtx.p <- edf.mtx.m
    }
    head(edf.mtx.p)
    
    # allocate a matrix
    non.seizure.mtx <- matrix(NaN, nrow = nn * f * window, ncol = ncol(edf.mtx.p)) 
    colnames(non.seizure.mtx) <- colnames(edf.mtx.m)
    
    nrows.edf <- nrow(edf.mtx.p)
    mm <- matrix(NaN, nrow = nn, ncol = 4)
    for (k in 1:nn) {
      if (random) {
        r <- sample(1:(nrows.edf - f * window), 1) # randomly 
      } else {
        set.seed(my.seed + k)
        r <- sample(1:(nrows.edf - f * window), 1) # seed
      }  
      temp <- as.matrix(edf.mtx.p[r:(r + f * window - 1), ])
      
      if (filtering) txt = "f" else txt = ""
      ext <- ".txt"
      file = paste(
        dir, 
        "working/bin_files/",
        txt,
        "ns", 
        "_e", which.expert, 
        "_p", formatC(q, width = 2, format = "d", flag = "0"),  
        "_w",  formatC(window, width = 2, format = "d", flag = "0"),
        "c", formatC(chunks, width = 2, format = "d", flag = "0"),
        "_seq_", formatC(k, width = 4, format = "d", flag = "0"),
        "_",
        f, "Hz", sep = "") 
      file.txt <- paste(file, ext, sep = "")
      
      # Write to a text file (just to make it easier to view the content, we don't really need this file)
      if (write.bin.as.txt.files) {
        write.table(
          x = round(temp, 4),
          file = file.txt,
          dec = ".",
          sep = "\t",
          col.names = TRUE,
          row.names = FALSE
        )
      }
      
      # To binary file
      # https://gregstoll.com/~gregstoll/floattohex/
      if (write.bin.files) {
        ext = ".bin"
        file.bin <- paste(file, ext, sep = "")
        fh = file(file.bin, "wb")
        for (n in 1:nrow(temp)) {
          # "-2", because the last two columns contain the timestamp and class label
          data_row <- as.numeric(temp[n, 1:(ncol(edf.mtx.p) - 2)])
          writeBin(data_row, fh) 
        }
        close(fh)
      }
      
      non.seizure.mtx[((k - 1) * f * window + 1):(k * f * window),] <- temp
      mm[k,] = t(c(q, k, r, (r + f * window - 1)))
    } # for (k in 1:nn)
    
    head(non.seizure.mtx)

    NON.SEIZURE <- rbind(NON.SEIZURE, non.seizure.mtx)
    
    
    colnames(mm) <- c("patient", "sequence_number", "from_sample",	"to_sample")
    
    if (q == non.seizure.IDs[1]) cn = TRUE else cn = FALSE
    
    # We are temporarily disabling warnings to avoid displaying an unobtrusive warning
    options(warn = -1)
    
    write.table(
      x = mm,
      file = paste(dir, "working/aux_files/non_seizures_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      append = TRUE,
      sep = "\t",
      col.names = cn,
      row.names = FALSE
    )
    
    options(warn = 0)
    
    colnames(NON.SEIZURE) <- colnames(edf.mtx.m)
    head(NON.SEIZURE)
    
  } # for (ii in non.seizure.IDs)
  
  # In the last column 0 - no seizure 
  NON.SEIZURE[, cc] <- rep(0, nrow(NON.SEIZURE))
  colnames(NON.SEIZURE)[cc] <- "seizure"
  head(NON.SEIZURE)
  nrow(NON.SEIZURE)
  
  if (write.txt.files) {
    write.table(
      x = round(NON.SEIZURE, 4),
      file = paste(dir, "working/aux_files/NON.SEIZURE_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""), 
      dec = ".",
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  }
  
  # ///////////////////////////////////////////////////////////////////////////////////////////
  # Combining the partial results into the final one file 
  # ///////////////////////////////////////////////////////////////////////////////////////////
  FINAL <- rbind(SEIZURE, NON.SEIZURE)
  FINAL.mtx <- as.matrix(FINAL)
  
  if (write.hdf5.files) {
    fname <- paste(dir, "working/aux_files/expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.hdf5", sep ="")
    unlink(fname)
    h5write(obj = FINAL.mtx, file = fname, name = "FINAL.mtx", createnewfile = TRUE)
    h5closeAll()
  }  

  if (write.txt.files) {
    write.table(
      x = round(FINAL, 4),
      file = paste(dir, "working/aux_files/expert_", which.expert, "_", window, "sec_", chunks, "chunk_", f, "Hz.txt", sep = ""),
      sep = "\t", 
      append = FALSE,
      col.names = FALSE,
      row.names = FALSE
    )
  } 
  
  list(
    SEIZURE = SEIZURE,
    NON.SEIZURE =NON.SEIZURE,
    FINAL = FINAL
  )
}

# ///////////////////////////////////////////////////////////////////////////////////////////
mp2tf <- function(
  SQLiteFile, 
  channel,  
  mode = "sqrt",
  freqDivide = 1, 
  increaseFactor= 1,
  displayCrosses = TRUE,
  grid = FALSE,
  crossesColor = "white",
  palette = 'Lajolla', 
  rev = TRUE,
  outMode = "plot",
  fileName = NA,
  fileSize = c(512, 512),
  reportFile = "report",
  plotSignals = TRUE) {
  
  # SQLiteFile         - plik, który "schodzi" z programu empi
  # channel            - który kanał z pliku sqlite obrabiać
  # mode               - "sqrt", "log", "linear" (b. podobnie jak w Svarog)
  # freqDivide         - ile razy od góry ograniczyć wyświetlaną częstotliwość. Np. Gdy f=256Hz
  #                      to maksymalna częstotliwość jest f/2 (reguła Nyquista) a po ograniczeniu
  #                      mamy f/2/freqDivide
  # increaseFactor     - współczynnik zwiększenie ilości pikseli w osi f. Robimy tak, gdyż 
  #                      po ograniczeniu od góry częstotliwości parametrem freqDivide w blobach
  #                      zaczyna byc widać pikselizację, co nie wygląda ładnie.
  #                      Można podawać dowolne wartości ale najbardziej sensowne są liczby
  #                      całkowite dodanie (np. 2, 4, 5, 8)
  # displayCrosses     - czy w środkach atomów maja być wyświetlone małe krzyżyki
  # palette            - paleta z listy zwracanej przez hcl.pals() 
  #                      lub 'my custom palette' (dale kolory bardo podobne do Svarog-a)
  # rev                - rev param in hcl.colors
  # outMode            - 'plot': rysuje na ekranie, 
  #                      'file': zapisuje obraz na pliku 'fileName'
  #                      'RData': zapisuje dane potrzebne do utworzenia rysunku w rozmiarze fileSize
  #                      (w pliku 'fileName' z rozszerzenie .RData)
  #                      'console': tylko do pliku zapisuje podstawowe dane odczytane z plików db
  # reportFile         - nazwa pliku, w który zapisane będą dane o ilości atomów powstałych
  #                      w wyniku działania programu empi.exe
  # plotSignals        - czy mają być wyświetlane też sygnały 'original' oraz 'reconstructed'
  
  if (outMode != "plot" & outMode != "file" & outMode != "RData" & outMode != "console")
    stop("\n--> Incorrect value for 'outMode' parameter' <--")
  
  
  if (palette == 'my custom palette') {
    col <-  c(
      "#000f82", "#001385", "#011789", "#011b8d", "#021f91", "#022395", "#032798", "#042b9c",
      "#042fa0", "#0533a4", "#0537a8", "#063bab", "#073faf", "#0743b3", "#0847b7", "#084bbb",
      "#094fbf", "#0a53c2", "#0a57c6", "#0b5bca", "#0b5fce", "#0c63d2", "#0d67d5", "#0d6bd9",
      "#0e6fdd", "#0e73e1", "#0f77e5", "#107be8", "#107fec", "#1183f0", "#1187f4", "#128bf8",
      "#1390fc", "#1693f8", "#1996f4", "#1c9af0", "#1f9dec", "#22a1e8", "#25a4e4", "#28a8e0",
      "#2cabdc", "#2faed8", "#32b2d4", "#35b5d0", "#38b9cc", "#3bbcc8", "#3ec0c4", "#41c3c0",
      "#45c7bc", "#48cab8", "#4bcdb4", "#4ed1b0", "#51d4ac", "#54d8a8", "#57dba4", "#5adfa0",
      "#5ee29c", "#61e598", "#64e994", "#67ec90", "#6af08c", "#6df388", "#70f784", "#73fa80",
      "#77fe7c", "#7bf978", "#7ff575", "#83f171", "#88ed6e", "#8ce96a", "#90e567", "#94e063",
      "#99dc60", "#9dd85c", "#a1d459", "#a5d055", "#aacc52", "#aec74e", "#b2c34b", "#b6bf47",
      "#bbbb44", "#bfb741", "#c3b33d", "#c7af3a", "#ccaa36", "#d0a633", "#d4a22f", "#d89e2c",
      "#dd9a28", "#e19625", "#e59121", "#e98d1e", "#ee891a", "#f28517", "#f68113", "#fa7d10",
      "#ff790d", "#fb750c", "#f8720c", "#f56f0b", "#f26c0b", "#ef690a", "#eb660a", "#e8630a",
      "#e56009", "#e25c09", "#df5908", "#db5608", "#d85308", "#d55007", "#d24d07", "#cf4a06",
      "#cc4706", "#c84306", "#c54005", "#c23d05", "#bf3a04", "#bc3704", "#b83404", "#b53103",
      "#b22e03", "#af2a02", "#ac2702", "#a82402", "#a52101", "#a21e01", "#9f1b00", "#9c1800",
      "#991500")
  } else {
    col <- hcl.colors(128, palette, rev = rev)
  }
  
  
  lDataFrames <- readSQLite(SQLiteFile)
  
  # sampling rate in Hz
  f <- as.numeric(lDataFrames[[2]][["value"]][3]) 
  
  # number of samples
  epochSize <- lDataFrames[[4]][["sample_count"]] 
  s <- epochSize / f # ilość sekund
  
  # grid size in t 
  t <- seq(from = 0, to = s , by = 1 / f) 
  
  # according to the Nyquist criterion
  maxf <- round((f / 2) / freqDivide) 
  
  # grid size in f
  y <- seq(from = 0, to = maxf, by = 1 / increaseFactor)  
  
  # in the empi program channels are numbered from 0
  rows <- which(lDataFrames[[1]]$channel_id == (channel - 1)) 
  
  # correction factor
  dd <- epochSize / f
  
  # t-f map
  tf.map <- matrix(0, nrow = epochSize, ncol = maxf * increaseFactor)
  
  grid.x <- seq(from = 0, to = s, length.out = epochSize)
  grid.y <- seq(from = 0, to = f / 2 / freqDivide, length.out = maxf * increaseFactor)
  mtx <- as.matrix(expand.grid(grid.x, grid.y))
  
  if (length(rows) == 0) {
    stop("\n--> There is no channel number ", channel,  " <--", sep = "")
  }
  
  if (outMode == "console") {
    data <- paste(SQLiteFile, "\t", "channel: ", channel, "\t", "Number of atoms: ", length(rows), sep = "")
    write.table(data, file = reportFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    return(1)
  }
  
  cat("\n\nPlik SQLite: ", SQLiteFile, "\n", sep = "")
  cat("Channel: ", channel, "\n", sep = "")
  cat("Number of atoms: ", length(rows), "\n", sep = "")
  cat("Sampling rate: ", f, "\n", sep = "")
  cat("Epoch size (in points): ", epochSize, "\n", sep = "")
  cat("Signal length (in seconds): ", s, "\n", sep = "")

  data <- paste(SQLiteFile, "\t", "channel: ", channel, "\t", "number of atoms: ", length(rows), sep = "")
  write.table(data, file = reportFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  temp <- lDataFrames[[3]][["samples_float32"]][channel]
  utemp <- (unlist(temp))
  
  # number of atoms
  num.atoms <- length(rows)
  
  # parameters of individual atoms
  amplitude <- as.vector(as.matrix(lDataFrames[[1]][4][rows,]))
  energy <- as.vector(as.matrix(lDataFrames[[1]][5][rows,]))
  envelope <- as.vector(as.matrix(lDataFrames[[1]][6][rows,]))
  frequency <- as.vector(as.matrix(lDataFrames[[1]][7][rows,])) 
  phase <- as.vector(as.matrix(lDataFrames[[1]][8][rows,]))
  scale <- as.vector(as.matrix(lDataFrames[[1]][9][rows,])) 
  position <- as.vector(as.matrix(lDataFrames[[1]][10][rows,]))
  
  atoms <- list(
    amplitude = amplitude, 
    energy = energy, 
    envelope = envelope,
    frequency = frequency,
    phase = phase,
    scale = scale,
    position = position)
  
  if (outMode != "RData") {
    reconstruction <- rep(0, epochSize)
    gabor <- matrix(nrow = num.atoms, ncol = epochSize)
    
    for (i in 1:num.atoms) {
      out <- Gabor.fun(epochSize, f, position[i], phase[i], scale[i], frequency[i])
      reconstruction <- reconstruction + out$gabor * sqrt(energy[i] * f)
      gabor[i,] <- out$gabor
    }
    
    # We read the input data from the .db file (they are stored there as float32 numbers)
    # For example: c0 74 23 f3  =  -3.81469
    originalSignal <- NA
    
    for (i in 1:epochSize) {
      (b <- readBin(utemp[((i - 1) * 4 + 1) : ((i - 1) * 4 + 4)], "raw", 4))
      # swap to use big-endian
      (b2 <- paste(b[4], b[3], b[2], b[1], sep = ""))
      # https://stackoverflow.com/questions/39461349/converting-hex-format-to-float-numbers-in-r
      originalSignal[i] <- readBin(as.raw(strtoi(substring(b2, (step <- seq(1, nchar(b2), by = 2)), step + 1), 16)), "double", n = 1, size = 4)
    }
  } else {
    reconstruction <- NA
    originalSignal <- NA
    gabor <- NA
  }
  
  # https://en.wikipedia.org/wiki/Error_function
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1 
  SQRT_PI = sqrt(pi)
  BI_SQRT_PI = 2 * SQRT_PI
  tx <- seq(1:(length(t) - 1))
  fy <- seq(1:(length(y) - 1))
  tlowx <- tx - 1
  tuppx <- tx
  flowy <- (fy - 1) / increaseFactor
  fuppy <- fy / increaseFactor
  gauss.tx <- matrix(0, nrow = num.atoms, ncol = (length(t) - 1))
  gauss.fy <- matrix(0, nrow = num.atoms, ncol = maxf * increaseFactor)

  for (n in 1:num.atoms) {

    gt1 <- (BI_SQRT_PI / (scale[n] * f)) * (tlowx - (position[n] * f))
    gt2 <- (BI_SQRT_PI / (scale[n] * f)) * (tuppx - (position[n] * f))
    gauss.tx[n,] <- (erf(gt2) - erf(gt1)) 

    gf1 <- (SQRT_PI * (scale[n] * f) / epochSize) * dd * (flowy  - (frequency[n]))
    gf2 <- (SQRT_PI * (scale[n] * f) / epochSize) * dd * (fuppy  - (frequency[n]))
    gauss.fy[n,] <- (erf(gf2) - erf(gf1)) 

    if (mode == "sqrt") 
      tf.map <- tf.map + kronecker(gauss.tx[n,], t(gauss.fy[n,])) * sqrt(energy[n] * f)
    
    if (mode == "log") 
      tf.map <- tf.map + kronecker(gauss.tx[n,], t(gauss.fy[n,])) * log(energy[n] * f)
    
    if (mode == "linear") 
      tf.map <- tf.map + kronecker(gauss.tx[n,], t(gauss.fy[n,])) * energy[n] * f
      
  } # for (n in 1:num.atoms)
  
  if (outMode == "plot") {  
    if (plotSignals) {
      grid.matrix <- cbind(c(1,1,1,2,3))
      layout(grid.matrix, widths = c(2,2,2), heights = c(3,1,1))
      # mai: c(bottom, left, top, right)
      par(pty = "m", mai = c(0.55, 0.6, 0.2, 0.4))
    } else {
      par(mfrow = c(1,1), pty = "m")
      par(mai = c(0.9, 0.9, 0.2, 0.4))
    }
    
    graphics::image(
      x = t, y = y, z = tf.map,
      xlab = "Time [s]", ylab = "Frequency [Hz]", cex.axis = 1.0, cex.lab = 1.0, col = col
    )

    # W środkach atomów wyświetlamy małe krzyżyki
    if (displayCrosses) 
      points(position, frequency, pch = 3, col = crossesColor, cex = 0.7)
    
    if (grid) 
      grid(col = "grey")
    
    if (plotSignals) {
      xx <- seq(from = 0, to = epochSize / f, length.out =  epochSize)
      plot(x = xx, originalSignal, type = "l", xlab = "", ylab = "", xaxs = "i", main = "original", panel.first = grid())
      abline(h = 0, col = "blue")
      
      plot(x = xx, reconstruction, type = "l", xlab = "", ylab = "", xaxs = "i", main = "reconstructed", panel.first = grid())
      abline(h = 0, col = "blue")
    }
    
  } # if (outMode == "plot")
  
  if (outMode == "file") {
    graphics.off()
    png(fileName, width = fileSize[1], height = fileSize[2], pointsize = 18)
    par(pty = "m", mai = c(0, 0, 0, 0))
    
    
    graphics::image(
      x = t, y = y, z = tf.map,
      xlab = "", ylab = "", col = col,
      yaxt = "n", xaxt = "n"
    )
    
    if (displayCrosses) 
      points(position, frequency, pch = 3, col = crossesColor, cex = 0.7)      
    
    dev.off()
  } # if (outMode == "file")
  
  if (outMode == "RData") { 
    zz  <- tf.map      
    # Usuwamy rozszerzenie
    fileName <- paste(tools::file_path_sans_ext(fileName), ".RData", sep = "")
    
    rr <- raster::raster(nrow = ncol(zz), ncol = nrow(zz)) # tu musi właśnie być odwrotnie
    rr[] <- t(zz)
    tt <- raster::raster(ncol = fileSize[1], nrow = fileSize[2])
    tt <- raster::resample(rr, tt)
    m2 <- matrix(tt@data@values, fileSize[1], fileSize[2])
    #image(m2, col = col)
    #range(m2)
    # Przeskalowanie do zakresu 0-1
    # Zabezpieczamy się, przed sytuacją, gdy w mianowniku pojawia się zero
    if (max(m2) - min(m2) == 0) {
      tf.matrix <- matrix(0, fileSize[1], fileSize[2])  
    } else {
      tf.matrix <- (m2 - min(m2)) / (max(m2) - min(m2))
    }
    save(tf.matrix, file = fileName)
  } # if (outMode == "RData")
  
  list(
    gaborFunctions = gabor, 
    atoms = atoms,
    reconstructedSignal = reconstruction,
    originalSignal = originalSignal,
    f = f,
    epochSize = epochSize,
    numberOfSecs = s,
    z = tf.map,
    freqDivide = freqDivide)
}

# ///////////////////////////////////////////////////////////////////////////////////////////
readSQLite <- function(SQLite_db_file) {
  con <- dbConnect(drv = RSQLite::SQLite(), dbname = SQLite_db_file)
  
  ## list all tables
  tables <- dbListTables(con)
  
  ## create a data.frame for each table
  lDataFrames <- vector("list", length = length(tables))
  for (i in seq(along = tables)) {
    lDataFrames[[i]] <- dbGetQuery(conn = con, statement = paste("SELECT * FROM '", tables[[i]], "'", sep = ""))
  }
  
  dbDisconnect(con)
  return(lDataFrames)
}

# ///////////////////////////////////////////////////////////////////////////////////////////
Gabor.fun <- function(
  numberOfSamples, 
  samplingFrequency, 
  mean, 
  phase, 
  sigma, 
  frequency, 
  normalization = T) {
  
  omega <- 2 * pi * frequency
  t <- seq(from = 0, to = numberOfSamples - 1, by = 1) / samplingFrequency
  v <- exp(-pi * ((t - mean) / sigma)^2)
  u <- cos(omega * (t - mean) + phase)
  gabor <- u * v
  if (normalization) {
    gabor <- vec.norm(gabor)
  }
  list(cosinus = u, gauss = v, gabor = gabor, t = t)
}

# ///////////////////////////////////////////////////////////////////////////////////////////
vec.norm <- function(x) {x / sqrt(sum(x^2))}

# ///////////////////////////////////////////////////////////////////////////////////////////
generate_RData_files <- function(
  pattern, 
  tf_maps_destination, 
  tf_file_size, 
  sqlite_source, 
  report_file, 
  progress_file, 
  fs, 
  fns) {
  
  files <- list.files(sqlite_source, include.dirs = TRUE, pattern = pattern)
  # Remove extension
  files <- tools::file_path_sans_ext(files)
  
  fs <- fs
  fns <- fns
  
  if (pattern == 'fs') k <- fs
  if (pattern == 'fns') k <- fns
  
  for (i in k:length(files)) {
    
    g <- grep('fs', files[i])
    if (length(g) == 1) {
      format.str <- formatC(fs, width = 5, format = "d", flag = "0")
      fs <- fs + 1
    } 
    if (length(g) == 0) {
      format.str <- formatC(fns, width = 5, format = "d", flag = "0")
      fns <- fns + 1
    } 
    
    for (j in 1: 18) { # we have 18 channels in our files
      
      tfFileName <- paste(
        tf_maps_destination,
        files[i],
        "_ch",
        formatC(j, width = 2, format = "d", flag = "0"),
        "_",
        tf_file_size[1],
        "_",
        tf_file_size[2],
        "_file",
        format.str,
        ".png", sep = "")
      
      out <- mp2tf(
        SQLiteFile = paste(sqlite_source, files[i], ".db", sep = ""),
        channel = j,
        mode = "sqrt",
        # 4 means we are limited to 32Hz
        freqDivide = 4,
        increaseFactor= 16,
        displayCrosses = FALSE,
        grid = FALSE,
        crossesColor = "white",
        outMode = "RData",
        fileName = tfFileName,
        fileSize = tf_file_size,
        reportFile =report_file
      )
    }
    df <- data.frame(fs = fs, fns = fns, i = i, total =  length(files), row.names = NULL)
    write.table(df, progress_file, row.names = F, col.name = TRUE, sep = "\t", quote = FALSE)    
  }
}  

# ///////////////////////////////////////////////////////////////////////////////////////////
signal_energy <- function(SQLiteFile, channel, verbose = FALSE) {
  lDataFrames <- readSQLite(SQLiteFile)
  
  # Sampling rate in Hz
  f <- as.numeric(lDataFrames[[2]][["value"]][3]) 
  
  # In empi we number channels from 0
  rows <- which(lDataFrames[[1]]$channel_id == (channel - 1)) 
  
  # Number of atoms
  num.atoms <- length(rows)
  
  # Parameters of individual atoms
  amplitude <- as.vector(as.matrix(lDataFrames[[1]][4][rows,]))
  energy <- as.vector(as.matrix(lDataFrames[[1]][5][rows,]))
  envelope <- as.vector(as.matrix(lDataFrames[[1]][6][rows,]))
  frequency <- as.vector(as.matrix(lDataFrames[[1]][7][rows,])) 
  phase <- as.vector(as.matrix(lDataFrames[[1]][8][rows,]))
  scale <- as.vector(as.matrix(lDataFrames[[1]][9][rows,])) 
  position <- as.vector(as.matrix(lDataFrames[[1]][10][rows,]))
  
  # Number of samples
  epochSize <- lDataFrames[[4]][["sample_count"]] 
  
  if (verbose) {
    cat("\n\nPlik SQLite: ", SQLiteFile, "\n", sep = "")
    cat("Kanał: ", channel, "\n", sep = "")
    cat("Ilość atomów: ", num.atoms, "\n", sep = "")
  }
  
  # Based on the results returned by empi, we reconstruct the individual atoms
  reconstructedSignal <- rep(0, epochSize)
  gabor <- matrix(nrow = num.atoms, ncol = epochSize)
  
  for (i in 1:num.atoms) {
    out <- Gabor.fun(epochSize, f, position[i], phase[i], scale[i], frequency[i])
    reconstructedSignal <- reconstructedSignal + out$gabor * sqrt(energy[i] * f)
    gabor[i,] <- out$gabor
  }
  
  temp <- lDataFrames[[3]][["samples_float32"]][channel]
  utemp <- (unlist(temp))
  
  # We read the input data from the .db file (they are stored there as float32 numbers)
  # For example: c0 74 23 f3  =  -3.81469
  originalSignal <- NA
  
  for (i in 1:epochSize) {
    (b <- readBin(utemp[((i - 1) * 4 + 1) : ((i - 1) * 4 + 4)], "raw", 4))
    # swap to use big-endian
    (b2 <- paste(b[4], b[3], b[2], b[1], sep = ""))
    # https://stackoverflow.com/questions/39461349/converting-hex-format-to-float-numbers-in-r
    originalSignal[i] <- readBin(as.raw(strtoi(substring(b2, (step <- seq(1, nchar(b2), by = 2)), step + 1), 16)), "double", n = 1, size = 4)
  }
  
  # Signal energy
  o <- round(sum(originalSignal^2), 2)
  r <- round(sum(reconstructedSignal^2), 2)
  
  if (verbose) {
    cat("Energy of the original signal: ",o , "\n", sep = "")
    cat("Signal energy after reconstruction: ",r , "\n", sep = "")
    cat("reconstruction / original %: ", r / o * 100, "\n", sep = "")
  }
  
  list(
    SQLiteFile = SQLiteFile, 
    channel = channel,
    energy.OriginalSignal = o,
    energy.ReconstructedSignal = r,
    r.o.ratio = r / o * 100
  )
}    






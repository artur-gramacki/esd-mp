# file __run_it_first__.R
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

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --- IMPORTANT NOTE 1 ---
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# All codes were run and tested in version 4.4.1.
# When running scripts in version 4.5.1, a sudden system crash was observed. 
# We were unable to find the cause.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --- IMPORTANT NOTE 2 ---
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Before starting, download the dataset of neonatal EEG recordings to the 'edf' directory.  
# These are available at https://zenodo.org/record/4940267. 
# There are 79 EDF files and 3 CSV annotations files. 
# The EDF files are approximately 4GB in size.  

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --- IMPORTANT NOTE 3 ---
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Before starting, make sure that the 'R' directory is the current directory,
# as returned from the getwd() function.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --- IMPORTANT NOTE 4 ---
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Set the correct path to pyton.exe executable file

dir <- "../"
python_dir <- "c:/Programy/miniconda3/"

# Set to TRUE to remove the given directories (if they are present) and recreate them
recreate_dirs <- TRUE

if(recreate_dirs) {
  unlink(paste(dir, 'working/', sep = ""), recursive = TRUE, force = TRUE)
  dir.create(paste(dir, 'working/', sep = ""))
  dir.create(paste(dir, 'working/aux_files', sep = ""))
  dir.create(paste(dir, 'working/bin_files', sep = ""))
  dir.create(paste(dir, 'working/bin_files/eA_w10_c20', sep = ""))
  dir.create(paste(dir, 'working/bin_files/eB_w10_c20', sep = ""))
  dir.create(paste(dir, 'working/bin_files/eC_w10_c20', sep = ""))
  dir.create(paste(dir, 'working/hdf5_files', sep = ""))
  dir.create(paste(dir, 'working/hdf5_files_STFT', sep = ""))
  dir.create(paste(dir, 'working/sqlitedb_files', sep = ""))
  dir.create(paste(dir, 'working/sqlitedb_files/eA_w10_c20', sep = ""))
  dir.create(paste(dir, 'working/sqlitedb_files/eB_w10_c20', sep = ""))
  dir.create(paste(dir, 'working/sqlitedb_files/eC_w10_c20', sep = ""))
  dir.create(paste(dir, 'working/tf_maps', sep = ""))
  dir.create(paste(dir, 'working/tf_maps/eA_w10_c20_64_64', sep = ""))
  dir.create(paste(dir, 'working/tf_maps/eB_w10_c20_64_64', sep = ""))
  dir.create(paste(dir, 'working/tf_maps/eC_w10_c20_64_64', sep = ""))
  dir.create(paste(dir, 'working/tf_maps_STFT', sep = ""))
  dir.create(paste(dir, 'working/tf_maps_STFT/eA_w10_c20_64_64', sep = ""))
  dir.create(paste(dir, 'working/tf_maps_STFT/eB_w10_c20_64_64', sep = ""))
  dir.create(paste(dir, 'working/tf_maps_STFT/eC_w10_c20_64_64', sep = ""))
}

# After completing this script:
# 1. working/aux_files directory should contain 9 files
# 2. working/bin_files/eA_w10_c20 directory should contain 48 files
# 3. working/bin_files/eB_w10_c20 directory should contain 82 files
# 4. working/bin_files/eC_w10_c20 directory should contain 66 files
source("edf_to_bin.R")

# After completing this script:
# 1. working/sqlitedb_files/eA_w10_c20 directory should contain 48 files
# 2. working/sqlitedb_files/eB_w10_c20 directory should contain 82 files
# 3. working/sqlitedb_files/eC_w10_c20 directory should contain 66 files
source("bin_to_db.R")

# After completing this script:
# 1. working/tf_maps/eA_w10_c20_64_64 directory should contain 864 files
# 2. working/tf_maps/eB_w10_c20_64_64 directory should contain 1476 files
# 3. working/tf_maps/eC_w10_c20_64_64 directory should contain 1188 files
# 4. working/tf_maps directory should contain 12 log/diagnostic files
source("db_to_RData.R")

# After completing this script:
# 1. working/tf_maps_STFT/eA_w10_c20_64_64 directory should contain 864 files
# 2. working/tf_maps_STFT/eB_w10_c20_64_64 directory should contain 1476 files
# 3. working/tf_maps_STFT/eC_w10_c20_64_64 directory should contain 1188 files
source("bin_to_RData_STFT.R")

# After completing the below command:
# 1. working/hdf5_files directory should contain 4 files:
#   eABC_w10_c20_64_64_XY.hdf5
#   eA_w10_c20_64_64.hdf5
#   eB_w10_c20_64_64.hdf5
#   eC_w10_c20_64_64.hdf5
# 2. working/hdf5_files_STFT directory should contain 4 files:
#   eABC_w10_c20_64_64_XY.hdf5
#   eA_w10_c20_64_64.hdf5
#   eB_w10_c20_64_64.hdf5
#   eC_w10_c20_64_64.hdf5
run_command  <- paste(python_dir, "python.exe ", " ../Python/RData_to_hdf5.py", sep = "")
system(run_command)

# After completing the below command:
# 1. working/hdf5_files directory should contain additional 60 files (5 folds for train, validate, test):
# 2. working/hdf5_files_STFT directory should contain additional 60 files (5 folds for train, validate, test):
run_command  <- paste(python_dir, "python.exe ", " ../Python/split_data.py", sep = "")
system(run_command)

############################################################################################
############################################################################################
# Demonstrations. Some examples showing the most important elements of data processing
############################################################################################
############################################################################################

source("my_custom_palette.R")
source("funs.R")

############################################################################################
# MP decomposition of a sample EEG data (saved as 64-bit binary file)
# Sampling frequency was 256HZ, the data is 10s long
# So, the size of the sample.bin file is: 256 x 10 x 8 = 20480 bytes
############################################################################################
command <- paste(
  "../empi_1.0.3/empi.exe ",
  "sample/sample.bin ",
  "sample/sample.db ",
  "-f 256 -c 1 --channels 1 -o local --gabor -i 50 --input64",
  sep = "")
system(command)

############################################################################################
# Display sample time-frequency map on the screen
############################################################################################
out <- mp2tf(
  SQLiteFile = "sample/sample.db", 
  channel = 1, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 4,
  displayCrosses = TRUE, 
  grid = FALSE, 
  crossesColor = "white", 
  palette = "my custom palette", 
  rev = TRUE,
  outMode = "plot"
)

############################################################################################
# Save sample time-frequency map to the png file
############################################################################################
out <- mp2tf(
  SQLiteFile = "sample/sample.db", 
  channel = 1, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 16,
  displayCrosses = TRUE, 
  crossesColor = "white", 
  grid = FALSE, 
  palette = "my custom palette", 
  outMode = "file",
  fileName = "sample/sample.png",
  fileSize = c(512, 512),
)

############################################################################################
# Save time-frequency map to RData file
############################################################################################
out <- mp2tf(
  SQLiteFile = "sample/sample.db",
  channel = 1, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 16,
  displayCrosses = FALSE, 
  grid = FALSE, 
  palette = "my custom palette", 
  outMode = "RData",
  fileName = "sample/sample.RData",
  fileSize = c(64, 64),
)

############################################################################################
# Read RData file form disk and display it's content on the screen
# t-f map is 64x64 pixels in size. Then save it as png file. 
############################################################################################
load("sample/sample.RData")
par(pty = "s")
graphics::image(
  x = seq(0, 10, length.out = 64), 
  y = seq(0, 32, length.out = 64), 
  z = tf.matrix, 
  col = my_custom_palette,
  xlab = "Time [s]",
  ylab = "Frequency [Hz]" )

graphics.off()
png("sample/sample_64x64pixels.png", width = 64, height = 64, units = "px", type = "cairo-png")
par(pty = "m", mai = c(0, 0, 0, 0))
graphics::image(z = tf.matrix, col = my_custom_palette)
dev.off()

############################################################################################
# Read sample bin file, then execute STFT 
############################################################################################
fn <- "sample/sample.bin"
fs <- file.size(fn)
fh <- file(fn, "rb")
data <- NA
for (k in 1:(fs / 8)) { # 8: size of double in bytes
  b <- readBin(con = fh, what = "raw", n = 8)
  # swap to use little-endian
  hex_str <- paste(b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], sep = "")
  # https://stackoverflow.com/questions/39461349/converting-hex-format-to-float-numbers-in-r
  # Convert to raw
  raw_vec <- as.raw(strtoi(substring(hex_str, (step <- seq(1, nchar(hex_str), by = 2)), step + 1), 16))
  # Convert raw to double
  num <- readBin(raw_vec, what = "double", n = 1, size = 8, , endian = "little")
  data[k] <- num
}
close(fh)
# Much slower but gives t-f maps with better resolution in f axis
ft <- e1071::stft(X = data, win = 128, inc = 1, coef = 512)
#ft <- e1071::stft(X = data, win = 128, inc = 1, coef = 64)
# y axis limit to fmax Hz, our EEG recordings' sampling freq = 256
# 256 / 2: Nyquist freq
fmax <- 32
p <- (ncol(ft$values) * fmax) / (256 / 2)
par(pty = "s")
graphics::image(
  z = ft$values[,1:p], 
  col = my_custom_palette,
  xlab = "Time [s]",
  ylab = "Frequency [Hz]",
  xaxt = "n",
  yaxt = "n",
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

############################################################################################
# Create RData file from STFT results (64 x 64 pixels in size)
############################################################################################
zz <- ft$values[,1:p]      
rr <- raster::raster(nrow = ncol(zz), ncol = nrow(zz)) 
rr[] <- t(zz)
tt <- raster::raster(ncol = 64, nrow = 64)
tt <- raster::resample(rr, tt)
m2 <- matrix(tt@data@values, 64, 64)
# Rescaling to the range 0-1
# Protect against a situation where a zero appears in the denominator
if (max(m2) - min(m2) == 0) {
  tf.matrix <- matrix(0, 64, 64)  
} else {
  tf.matrix <- (m2 - min(m2)) / (max(m2) - min(m2))
}
save(tf.matrix, file = "sample/sample_STFT.RData")

############################################################################################
# Read Rdata file form disk and display it's content on the screen
# t-f map is 64x64 pixels in size. Then save it as png file. 
############################################################################################
load("sample/sample_STFT.RData")
par(pty = "s")
graphics::image(
  x = seq(0, 10, length.out = 64), 
  y = seq(0, 32, length.out = 64), 
  z = tf.matrix, 
  col = my_custom_palette,
  xlab = "Time [s]",
  ylab = "Frequency [Hz]" )

graphics.off()
png("sample/sample_STFT_64x64pixels.png", width = 64, height = 64, units = "px", type = "cairo-png")
par(pty = "m", mai = c(0, 0, 0, 0))
graphics::image(z = tf.matrix, col = my_custom_palette)
dev.off()

# file __RUN__.R
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
# Before starting, download the dataset of neonatal EEG recordings to the 'edf' directory.  
# These are available at https://zenodo.org/record/4940267. 
# There are 79 EDF files and 3 CSV annotations files. 
# The EDF files are approximately 4GB in size.  

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# --- IMPORTANT NOTE 2 ---
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Before starting, make sure that the 'edf_to_hdf5/R' directory is the current direcory,
# as returned from the getwd() function.

# ///////////////////////////////////////////////////////////////////////////////////////////
# Uncomment if you are sure you want to remove the given directories
# ///////////////////////////////////////////////////////////////////////////////////////////

dir <- "../"

# unlink(paste(dir, 'working/', sep = ""), recursive = TRUE)
# dir.create(paste(dir, 'working/', sep = ""))
# dir.create(paste(dir, 'working/aux_files', sep = ""))
# dir.create(paste(dir, 'working/bin_files', sep = ""))
# dir.create(paste(dir, 'working/bin_files/eA_w10_c20', sep = ""))
# dir.create(paste(dir, 'working/bin_files/eB_w10_c20', sep = ""))
# dir.create(paste(dir, 'working/bin_files/eC_w10_c20', sep = ""))
# dir.create(paste(dir, 'working/hdf5_files', sep = ""))
# dir.create(paste(dir, 'working/sqlitedb_files', sep = ""))
# dir.create(paste(dir, 'working/sqlitedb_files/eA_w10_c20', sep = ""))
# dir.create(paste(dir, 'working/sqlitedb_files/eB_w10_c20', sep = ""))
# dir.create(paste(dir, 'working/sqlitedb_files/eC_w10_c20', sep = ""))
# dir.create(paste(dir, 'working/tf_maps', sep = ""))
# dir.create(paste(dir, 'working/tf_maps/eA_w10_c20_64_64', sep = ""))
# dir.create(paste(dir, 'working/tf_maps/eB_w10_c20_64_64', sep = ""))
# dir.create(paste(dir, 'working/tf_maps/eC_w10_c20_64_64', sep = ""))

source("edf_to_bin.R")
source("bin_to_db.R")
source("db_to_RData.R")

# ///////////////////////////////////////////////////////////////////////////////////////////
# Some examples

source("my_custom_palette.R")

# Display time-frequency map on the screen
out <- mp2tf(
  SQLiteFile = "example.db", 
  channel = 1, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 16,
  displayCrosses = TRUE, 
  grid = FALSE, 
  crossesColor = "white", 
  palette = "my custom palette", 
  rev = TRUE,
  outMode = "plot"
)

# Save time-frequency map to file
out <- mp2tf(
  SQLiteFile = "example.db", 
  channel = 1, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 16,
  displayCrosses = TRUE, 
  crossesColor = "white", 
  grid = FALSE, 
  palette = "my custom palette", 
  outMode = "file",
  fileName = "example.png",
  fileSize = c(512, 512),
)

# Save time-frequency map to RData file
out <- mp2tf(
  SQLiteFile = "example.db", 
  channel = 1, 
  mode = "sqrt", 
  freqDivide = 4,
  increaseFactor = 16,
  displayCrosses = FALSE, 
  grid = FALSE, 
  palette = "my custom palette", 
  outMode = "RData",
  fileName = "example.RData",
  fileSize = c(64, 64),
)

# Read Rdata file form disk and display it's content on the screen
# t-f maps are 64x64 pixels size
load("example.RData")
par(pty = "s")
graphics::image(
  x = seq(0, 10, length.out = 64), 
  y = seq(0, 32, length.out = 64), 
  z = tf.matrix, 
  col = my_custom_palette,
  xlab = "t [s]",
  ylab = "f [Hz]" )

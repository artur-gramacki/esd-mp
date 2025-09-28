# file edf_to_bin.R
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

dir = "../"

# Check if edf files are already downloaded into the 'edf' dir.
files <- list.files(paste(dir, "edf", sep = ""), include.dirs = FALSE, pattern <- ".edf")
if (length(files) == 0) {
  stop(
  "Download to `edf` directory 79 EDF files of neonatal EEG recordings from https://zenodo.org/record/4940267. 
  See also `README.md` file.
  "
)
}

# Symbols of human experts
we <- c( "A", "B", "C")

# Annotations file names, as downloaded from https://zenodo.org/record/4940267
ann.f <- c("annotations_2017_A_fixed.csv",
          "annotations_2017_B.csv",
          "annotations_2017_C.csv")

# infant IDs which have seizures(s: seizure)
s.IDs <- c(1,4,5,7,9,11,13,14,15,16,17,19,20,21,22,25,31,34,36,38,39,40,41,44,47,50,51,52,62,63,66,67,69,71,73,75,76,77,78,79)
# for testitng it's better to use only a subset
s.IDs <- c(4)

# infant IDs which are seizure free (ns: not seizure)
ns.IDs <- c(3,10,18,27,28,29,30,32,35,37,42,45,48,49,53,55,57,58,59,60,70,72)
# for testitng it's better to use only a subset
ns.IDs <- c(10)

# To ensure repeatable results
my.seeds = c(42, 1024, 9999)

timestamp()
for (i in 1:3) {
  out <-  generate_samples(
    which.expert = we[i], 
    annotations_file = ann.f[i], 
    seizure.IDs = s.IDs, 
    non.seizure.IDs = ns.IDs, 
    window = 10, 
    chunks = 20, 
    down.sampling.factor = 1, 
    filtering = TRUE, 
    dir = dir, 
    random = FALSE, 
    my.seed = my.seeds[i],
    write.txt.files = FALSE, 
    write.hdf5.files = FALSE,
    write.bin.files = TRUE,
    write.bin.as.txt.files = FALSE
    )
  }
timestamp()

files <- list.files("../working/bin_files", include.dirs = TRUE, pattern = ".*eA.*.bin")
from_pathes <- file.path("../working/bin_files", files)
to_pathes <- file.path("../working/bin_files/eA_w10_c20",files)
out <- file.rename(from = from_pathes, to = to_pathes)

files <- list.files("../working/bin_files", include.dirs = TRUE, pattern = ".*eB.*.bin")
from_pathes <- file.path("../working/bin_files", files)
to_pathes <- file.path("../working/bin_files/eB_w10_c20",files)
out <- file.rename(from = from_pathes, to = to_pathes)

files <- list.files("../working/bin_files", include.dirs = TRUE, pattern = ".*eC.*.bin")
from_pathes <- file.path("../working/bin_files", files)
to_pathes <- file.path("../working/bin_files/eC_w10_c20",files)
out <- file.rename(from = from_pathes, to = to_pathes)

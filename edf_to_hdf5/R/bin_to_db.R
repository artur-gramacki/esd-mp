# file bin_to_db.R
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

patterns <- c(".*eA.*.bin", ".*eB.*.bin", ".*eC.*.bin")
dirs <- c("eA_w10_c20", "eB_w10_c20", "eC_w10_c20")  

for (i in 1:3) {
  files <- list.files(paste("../working/bin_files/", dirs[i], sep = ""), include.dirs = FALSE, pattern <- patterns[i])
  # remove the extension
  files <- tools::file_path_sans_ext(files)

  # From README.md:
  # -o local` enables local parameter optimization, but starting only from each iteration’s best match
  # and therefore, does not guarantee choosing the globally best atom in each iteration. However, it 
  # could constitute a sensible tradeoff between precision and performance for real-world usages, 
  # as it can be much faster than `-o global
  #              ^^^^^^^^^^^
  # "-r" and "--energy-error" get default values  (i.e. 0.001 and 0.05 respectively)
  
  for (j in 1:length(files)) {
    command <- paste(
      "../empi_1.0.3/empi.exe ",
      "../working/bin_files/",
      dirs[i],
      "/",
      files[j],
      ".bin ",
      "../working/sqlitedb_files/",
      dirs[i],
      "/",
      files[j],
      ".db ",
      "-f 256 -c 18 --channels 1-18 -o local --gabor -i 50 --input64",
      sep = "")
    system(command)
  }
}


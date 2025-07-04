# db_to_RData.R
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

dir <- "../"

generate_RData_files(
  pattern = "fns", 
  tf_maps_destination = paste(dir, "working/tf_maps/eA_w10_c20_64_64/", sep = ""), 
  tf_file_size = c(64, 64), 
  sqlite_source =  paste(dir, "working/sqlitedb_files/eA_w10_c20/", sep = ""),
  report_file = paste(dir, "working/tf_maps/report_A_fns.txt", sep = ""),
  progress_file = paste(dir, "working/tf_maps/progress_A_fns.txt", sep = ""),
  fs = 1,
  fns = 1)

generate_RData_files(
  pattern = "fs", 
  tf_maps_destination = paste(dir, "working/tf_maps/eA_w10_c20_64_64/", sep = ""), 
  tf_file_size = c(64, 64), 
  sqlite_source =  paste(dir, "working/sqlitedb_files/eA_w10_c20/", sep = ""),
  report_file = paste(dir, "working/tf_maps/report_A_fs.txt", sep = ""),
  progress_file = paste(dir, "working/tf_maps/progress_A_fs.txt", sep = ""),
  fs = 1,
  fns = 1)


generate_RData_files(
  pattern = "fns", 
  tf_maps_destination = paste(dir, "working/tf_maps/eB_w10_c20_64_64/", sep = ""), 
  tf_file_size = c(64, 64), 
  sqlite_source =  paste(dir, "working/sqlitedb_files/eB_w10_c20/", sep = ""),
  report_file = paste(dir, "working/tf_maps/report_B_fns.txt", sep = ""),
  progress_file = paste(dir, "working/tf_maps/progress_B_fns.txt", sep = ""),  
  fs = 1,
  fns = 1)

generate_RData_files(
  pattern = "fs", 
  tf_maps_destination = paste(dir, "working/tf_maps/eB_w10_c20_64_64/", sep = ""), 
  tf_file_size = c(64, 64), 
  sqlite_source =  paste(dir, "working/sqlitedb_files/eB_w10_c20/", sep = ""),
  report_file = paste(dir, "working/tf_maps/report_B_fs.txt", sep = ""),
  progress_file = paste(dir, "working/tf_maps/progress_B_fs.txt", sep = ""),
  fs = 1,
  fns = 1)

generate_RData_files(
  pattern = "fns", 
  tf_maps_destination = paste(dir, "working/tf_maps/eC_w10_c20_64_64/", sep = ""), 
  tf_file_size = c(64, 64), 
  sqlite_source =  paste(dir, "working/sqlitedb_files/eC_w10_c20/", sep = ""),
  report_file = paste(dir, "working/tf_maps/report_C_fns.txt", sep = ""),
  progress_file = paste(dir, "working/tf_maps/progress_C_fns.txt", sep = ""),  fs = 1,
  fns = 1)

generate_RData_files(
  pattern = "fs", 
  tf_maps_destination = paste(dir, "working/tf_maps/eC_w10_c20_64_64/", sep = ""), 
  tf_file_size = c(64, 64), 
  sqlite_source =  paste(dir, "working/sqlitedb_files/eC_w10_c20/", sep = ""),
  report_file = paste(dir, "working/tf_maps/report_C_fs.txt", sep = ""),
  progress_file = paste(dir, "working/tf_maps/progress_C_fs.txt", sep = ""),
  fs = 1,
  fns = 1)




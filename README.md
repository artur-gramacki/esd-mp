# esd-mp

This is the Electronic Supplements for the paper:

> ## A deep learning epileptic seizure detection based on matching pursuit algorithm and its time-frequency graphical representation ##

by *Mateusz Kunik & Artur Gramacki*

e-mails:  m.kunik@issi.uz.zgora.pl, a.gramacki@ck.uz.zgora.pl

A complete repository consists of:
1. R and Python source files for generating input files for neural network (`hdf5` files). All the required scripts and other files are stored in the [edt_to_hdf5](https://github.com/artur-gramacki/esd-mp/tree/main/edf_to_hdf5) directory. You should start by running the [\_\_RUN\_\_.R](https://github.com/artur-gramacki/esd-mp/blob/main/edf_to_hdf5/R/__RUN__.R) script. The [edf_to_hdf5/R]{https://github.com/artur-gramacki/esd-mp/blob/main/edf_to_hdf5/R/} directory must be the current directory. <br><br> IMPORTANT NOTE: the entire procedure for generating `hdf5` files takes a very long time, especially the [db_to_Rdata.R](https://github.com/artur-gramacki/esd-mp/blob/main/edf_to_hdf5/R/db_to_Rdata.R) script can take up to 2-3 days to execute. In the future we plan to rewrite this script in a low-level language (C, C++, Java) to shorten this execution time. However, for testing, there is no need to generate `hdf5` files from scratch. The [data](https://github.com/artur-gramacki/esd-mp/tree/main/data) directory contains links to the files already generated.

3. Python scripts for deep learning tasks. All the required scripts are stored in the [scripts](https://github.com/artur-gramacki/esd-mp/tree/main/scripts) directory.


In case of problems, the authors declare the necessary help for potential researchers.

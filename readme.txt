****************************************************************************************
                                       POODLE-S
                  Copyright (C) Kana Shimizu, All rights reserved, 2009
****************************************************************************************


How to use POODLE-S

1. Unzip the downloaded file.
2. cd POODLES-v1.0/src/
3. make
4. Edit execPoodleS.sh, and change the shell parameter "INSTALLDIR".
5. Edit execPoodleS_SRC, and change the shell parameter "INSTALLDIR", "LOCATION_BLASTPGP", "LOCATION_MAKEMAT", "LOCATION_NR"
6. Copy single fasta files to any directory you want to store the files. These fasta files must have extentions "fst".
7. Prepare a file which includes filenames of the fasta files without extensions as follows:
Filename_1
Filename_2
Filename_3
...
8. Execute POODLE-S as follows:
./execPoodleS.sh prepared_list_filename directory_where_you_store_fasta_files
9. The results are written in ./results/prepared_filename.res

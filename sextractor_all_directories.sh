time \ls *nd2 | parallel -j 1 "basename {} .nd2" | parallel --verbose -j 14 "cd {}; bash /home/jim/code/ngs_project/sextractor_directory.sh ."

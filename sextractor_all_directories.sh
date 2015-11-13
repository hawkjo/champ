time \ls *nd2 | parallel -j 1 "basename {} .nd2" | parallel --verbose -j 14 "cd {}; bash ../../../src/sextractor_directory.sh ."

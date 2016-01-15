miseq_info=`\ls *_R1_* | xargs gunzip -c | head -1 | awk -F ':' '{print $1 ":" $2 ":" $3}'`
echo $miseq_info
\ls *_R1_* | xargs gunzip -c | grep "^$miseq_info" | awk '{print $1}' | sed 's/^@//' > ../read_names/all_read_names.txt

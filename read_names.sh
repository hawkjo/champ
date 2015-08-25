if [ $# -eq 0 ]; then
    echo 'Usage: read_names.sh <list_of_files>'
    exit 1
fi
miseq_info=`\ls $* | xargs gunzip -c | head -1 | awk -F ':' '{print $1 ":" $2 ":" $3}'`
echo $miseq_info
\ls $* | xargs gunzip -c | grep "^$miseq_info" | awk '{print $1}' | sed 's/^@//' > read_names.txt

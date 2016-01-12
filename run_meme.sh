if [ $# -ne 1 ]; then
    echo "Usage: run_meme.sh <fpath>"
    exit 1;
fi

f=$1
out_dir="meme_results_`echo $f | sed 's/.fa*//'`"

meme $f -dna -o $out_dir -nostatus -time 18000 -maxsize 60000 -mod anr -nmotifs 3 -minw 6 -maxw 50 -revcomp

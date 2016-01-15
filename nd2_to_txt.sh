time ls /var/ngstest/*nd2 | parallel -j 14 python nd2_to_txt.py

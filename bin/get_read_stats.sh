#!/bin/bash

#echo -e "$@"

read_stats=$(zcat $@ | sed -n '2~4p' | perl -ne 'chomp; $cum=$cum+length($_); if (eof()){print "$cum ".sprintf("%.0f", $cum/$.)."\n"}')
#cum_length=$(echo -e "$read_stats" | cut -d " " -f 1)
#avg_length=$(echo -e "$read_stats" | cut -d " " -f 2)

echo -e "$read_stats"

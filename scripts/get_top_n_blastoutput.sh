#!/bin/bash
fname=$1
ntop=$2


# this is simple and short command
sort -k1,1 -k4nr -k3n $fname | awk -v top="$ntop" 'BEGIN{query=""; counter=1}{if(query!=$1){counter=1; print; query=$1;} else if(query==$1 && counter<top){print $0; counter+=1}}'

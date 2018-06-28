#!/bin/bash 

# using awk put the id and % together
# sort from 1 to 2 column in reverse order, meaning high to low
# again use, awk to select only the top record for each id as the top record has the highest % identity

if [ "$1" == "ascending" ]; then
	awk '{print $1, $3, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' $1 | sort -k1,2  -n | awk 'BEGIN{key=""}{if($1 != key){print; key=$1}}'
else if [ "$1" == "descending" ]; then
	awk '{print $1, $3, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' $1 | sort -k1,2 -r -n | awk 'BEGIN{key=""}{if($1 != key){print; key=$1}}'
else
	awk '{print $1, $3, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12}' $1 | sort -k1,2 -n | awk 'BEGIN{key=""}{if($1 != key){print; key=$1}}'
fi

#! /bin/bash


while IFS= read -r line; do
  input=$(echo "$line" | awk '{ printf "%17s, %17s, %17s,", $1, $2, $3 }')
  result=$(echo "$line" | cct -I +proj=cart +ellps=WGS84 | awk '{ printf "%15s, %14s, %5s", $1, $2, $3 }')
  echo "$input $result"
done

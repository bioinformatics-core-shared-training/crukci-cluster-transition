# Select lines from the middle of a file.
# Usage: middle.sh filename -end_line -num_lines
echo extract lines $(($2 - $3 + 1)) to $2 from file $1
head -$2 $1 | tail -$3 > middle_$1

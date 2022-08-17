# Filter regions file to only include lines where mean coverage was above and below specified threshold
zcat depressus_SRR5767284.regions.bed.gz | awk -F "\t" '{ if(($5 >= 10) && ($5 <= 100)) { print } }' - > dep_sco.bed

# Modify thresholds file to include proportions in addition to simple base counts at each threshold
zcat purpuratus_SRR7211988.thresholds.bed.gz | awk 'BEGIN {FS="\t" ; OFS="\t" } NR>1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $5+$6+$7+$8+$9, $5 ? ($5/($3-$2+1)) : 0, $6 ? ($6/($3-$2+1)) : 0, $7 ? ($7/($3-$2+1)) : 0, $8 ? ($8/($3-$2+1)) : 0, $9 ? ($9/($3-$2+1)) : 0}' - | gzip > test.bed.gz


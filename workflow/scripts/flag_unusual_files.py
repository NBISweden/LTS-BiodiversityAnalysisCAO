# This script checks the number of reads
# before and after subsetting. The goal
# is to flag samples that may have 
# contamination (e.g. human, bacteria etc).
# For non-contaminated samples, the number
# of subsetted reads should be less than
# or equal to the number of mapped and
# filtered reads.

# Run as:
# python flag_unusual_files.py all_filtered_reads.fasta subsetted_reads.fasta out.txt 


import os
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]
outfile=sys.argv[3]

# Get the number of lines in the original file
with open(f1, 'r') as og:
    for count1, line in enumerate(og):
        pass

# Get the number of lines in the subsetted file
with open(f2, 'r') as subsetted:
	for count2, line2 in enumerate(subsetted):
		pass

# Write summary to file
out = open(outfile, "w")
out.write('Original read file: ' + str(os.path.basename(f1)) + '\n')
out.write('Subsetted read file: ' + str(os.path.basename(f2)) + '\n')
out.write('Number of reads before subsetting:' + str((count1 + 1)/2) + '\n')
out.write('Number of reads in subset:' +  str((count2 + 1)/2) + '\n')
# Flag potentially contamined files
if count1 > count2:
	out.write("This file was subsetted and not all reads were used in the Blastn search. It could be a sign of contamination")
else:
	out.write('All reads were used in the Blastn search')
out.close()




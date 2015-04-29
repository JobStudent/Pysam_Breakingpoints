# Pysam_Breakingpoints
First script playing around pysam and finding breakingpoints 

__author__ = 'Xabier'
import pysam
from itertools import ifilter
from collections import Counter
samfile = pysam.AlignmentFile("break_points_hg19_sorted.bam", "rb")

chr1_i = samfile.references.index("chr1")
c = Counter()
from collections import defaultdict

def is_clipped(cigar_tuple):                                           # this is a function for finding the "clipped" reads
    return cigar_tuple[0] in (4, 5)

d=Counter()
e=Counter()
for i in d:
    d[i]+=1
for key, value in sorted(d.iteritems(), key=lambda (k,v): (v,k)):
    print "%s: %s" % (key, value)

for align in ifilter(lambda a:a.reference_id==chr1_i, samfile):       # a loop for searching "clipped" reads in the first chromosome
    c[is_clipped(align.cigar[0]) or is_clipped(align.cigar[-1])] +=1  # This line is irrelevant
    if is_clipped(align.cigar[0]):                                    # If the read starts with a clipped part, store in d dictionary, where
        d[align.pos]+=1                                               #  the key is the POSITION in the chromosome and the value is the NUMBER of reads (with a clipped part in the starting point) at that position in the chromosome.
    if is_clipped(align.cigar[-1]):                                   # Same as  d dictionary but in this case, we are working with reads that ends with a clipped point.
        e[align.aend]+=1

print "Start:", d.most_common()                                      # print dictionaries sorted by the value ( to see the positions in the chromosome with most clipped reads (This is suspicious of been breakingpoints)
print "End:", e.most_common()

samfile.close()

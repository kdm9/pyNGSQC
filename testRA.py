
##MEMORY LEAK IN LARGE FILES
import pyNGSQC
in_name = "/home/kevin/workspace/pel_huge.txt"
ra = pyNGSQC.FastqRandomAccess(in_name)
print len(ra.record_positions)
#for iii in xrange(len(ra.record_positions)):
#    print iii, ra.get(iii)[0]

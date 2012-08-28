import pyNGSQC
in_name = "/var/ws/borevitz/pel_1.small.fastq"
ra = pyNGSQC.FastqRandomAccess(in_name)
print ra.record_positions
for iii in xrange(len(ra.record_positions)):
    print iii, ra.get(iii)[0]

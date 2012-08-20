import pyNGSQC as ngs
fqrdr = ngs.FastqReader(
    "/home/kevin/UniWork/BIOL3157/Assignments/2/in.fastq",
    compression=ngs.NO_COMPRESSION
    )
print fqrdr
for read in fqrdr:
    print read

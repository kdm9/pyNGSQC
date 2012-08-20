import pyNGSQC as ngs
in_file = "/home/kevin/UniWork/BIOL3157/Assignments/2/in.fastq"
out_file = "/home/kevin/UniWork/BIOL3157/Assignments/2/out.fastq"
def test_1():
fqrdr = ngs.FastqReader(
    "/home/kevin/UniWork/BIOL3157/Assignments/2/in.fastq"
    )
print fqrdr
for read in fqrdr:
    print read

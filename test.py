import pyNGSQC as ngs
import QualFilter as qfil
in_file = "/home/kevin/UniWork/BIOL3157/Assignments/2/in.fastq"
out_file = "/home/kevin/UniWork/BIOL3157/Assignments/2/out.fastq"


def test_1():
    fqrdr = ngs.FastqReader(in_file)
    print fqrdr
    count = 0
    for read in fqrdr:
        count += 1
    print count
test_1()
print "Test 1 passed"


def test_2():
    qf = qfil.QualFilter(in_file, out_file)
    print qf
    qf.run()
test_2()
print "Test 2 passed"


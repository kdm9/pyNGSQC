import pyNGSQC as ngs
import QualFilter as qfil
import QualTrimmer as qtrim
import HardTrimmer as htrim
import BarcodeSplitter as bcs
import time
import csv
import unittest
#prefix = "E:/UniWork/BIOL3157/Assignments/2/"
#prefix = "/home/kevin/UniWork/BIOL3157/Assignments/2/"
prefix = "/home/kevin/workspace/"
#prefix = "/var/ws/borevitz/"
in_file = prefix + "in.fastq"  # "big.fastq" "pel_huge.txt"
out_dir = prefix + "out/"
timer = time.time()


class pyNGSQCtester(unittest.TestCase):
    timer = time.time()

    def testFastqReader(self):
        fqrdr = ngs.FastqReader(in_file)
        count = 0
        for read in fqrdr:
            count += 1
        print count
        print "testFastqReader took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(1, 1)

    def testQualFilter(self):
        qf = qfil.QualFilter(in_file, out_dir + "qf.fastq", qual_threshold=20)
        retval = qf.run()
        print "testQualFilter took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testQualTrimmer(self):
        qt = qtrim.QualTrimmer(in_file, out_dir + "qt.fastq", qual_threshold=20)
        retval = qt.run()
        print "testQualTrimmer took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testBarcodeSplitter(self):
        bc = bcs.BarcodeSplitter(in_file, out_dir)
        bc.set_barcodes_from_file(prefix + "barcodes.csv", csv.excel)
        retval = bc.run()
        print "testBarcodeSplitter took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testHardTrimmer(self):
        ht = htrim.HardTrimmer(in_file, out_dir + "ht.fastq", length=30)
        retval = ht.run()
        print "testHardTrimmer took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(pyNGSQCtester))
    return suite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=1).run(suite())

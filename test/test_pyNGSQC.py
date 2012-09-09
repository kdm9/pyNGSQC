#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pyNGSQC as ngs
import QualFilter as qfil
import QualStats as qstat
import QualTrimmer as qtrim
import HardTrimmer as htrim
import BarcodeSplitter as bcs
import time
import csv
import unittest
prefix = "./test/data/"
in_file = prefix + "in.fastq"
out_dir = prefix + "out/"
timer = time.time()


class pyNGSQCtester(unittest.TestCase):
    timer = time.time()

    def testFastqReader(self):
        fqrdr = ngs.FastqReader(in_file)
        count = 0
        for read in fqrdr:
            count += 1
        print "Count is: %i" % count
        print "testFastqReader took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(1, 1)

    def testQualFilter(self):
        qf = qfil.QualFilter(
            in_file,
            out_dir + "qf.fastq",
            qual_threshold=20,
            qual_offset=33
            )
        retval = qf.run_paralell()
        print "testQualFilter took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testQualTrimmer(self):
        qt = qtrim.QualTrimmer(
            in_file,
            out_dir + "qt.fastq",
            qual_threshold=20,
            min_length=10,
            qual_offset=33
            )
        retval = qt.run_paralell()
        print "testQualTrimmer took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testBarcodeSplitter(self):
        bc = bcs.ParalellBarcodeSplitter(in_file, out_dir)
        bc.set_barcodes_from_file(prefix + "barcodes.csv", csv.excel)
        retval = bc.run_paralell()
        print "testBarcodeSplitter took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testHardTrimmer(self):
        ht = htrim.HardTrimmer(in_file, out_dir + "ht.fastq", length=30)
        retval = ht.run_paralell()
        print "testHardTrimmer took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)

    def testQualStats(self):
        qs = qstat.QualStats(in_file, qual_offset=33)
        retval = qs.run()
        print "testQualStats took %.3f sec" % (time.time() - self.timer)
        self.timer = time.time()
        self.assertEqual(retval, True)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(pyNGSQCtester))
    return suite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=3).run(suite())

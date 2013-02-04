#!/usr/bin/env python
#Copyright 2012 Kevin Murray
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

import pyngsqc as ngs
from pyngsqc import qualfilter as qfil
from pyngsqc import qualstats as qstat
from pyngsqc import qualtrimmer as qtrim
from pyngsqc import hardtrimmer as htrim
from pyngsqc import barcodesplitter as bcs
from pyngsqc import collapser as col
from pyngsqc import converter as conv

import time
import csv
import unittest
import os

prefix = "./test/data/"
in_file = prefix + "test.fastq"
out_dir = prefix + "out/"
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)  # empty dir is not included in git, so make it here
timer = time.time()


class Tester(unittest.TestCase):

    def testFastqReader(self):
        fqrdr = ngs.FastqReader(in_file)
        count = 0
        for read in fqrdr:
            count += 1
        self.assertEqual(count, 1000)

    def testQualFilterParallel(self):
        qf = qfil.QualFilter(
            in_file,
            out_dir + "qf.fastq",
            qual_threshold=20,
            qual_offset=33,
            pass_rate=0.9,
            max_Ns=-1,
            )
        num_read, num_written = qf.run_parallel()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 996)

    def testQualFilter(self):
        qf = qfil.QualFilter(
            in_file,
            out_dir + "qf.fastq",
            qual_threshold=20,
            qual_offset=33,
            pass_rate=0.9,
            max_Ns=-1,
            )
        num_read, num_written = qf.run()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 996)

    def testCollapser(self):
        co = col.Collapser(
            in_file,
            out_dir + "col.fastq",
            )
        num_read, num_written = co.run()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 992)

    def testQualTrimmerParallel(self):
        qt = qtrim.QualTrimmer(
            in_file,
            out_dir + "qt.fastq",
            qual_threshold=20,
            min_length=10,
            qual_offset=33
            )
        num_read, num_written = qt.run_parallel()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 1000)

    def testQualTrimmer(self):
        qt = qtrim.QualTrimmer(
            in_file,
            out_dir + "qt.fastq",
            qual_threshold=20,
            min_length=10,
            qual_offset=33
            )
        num_read, num_written = qt.run()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 1000)

    def testBarcodeSplitterParallel(self):
        bc = bcs.BarcodeSplitter(in_file, out_dir, prefix + "barcodes.csv")
        num_read, num_written, barcode_counts = bc.run_parallel()
        expected_barcode_counts = {
                'CACACTTGAATC': 332,
                'CACCGATGTATC': 6,
                'CACATCACGATC': 7,
                'CACCAGATCATC': 281,
                'CACTTAGGCATC': 7,
                'CACGATCAGATC': 270,
                'CACGCCAATATC': 35,
                'CACTGACCAATC': 35,
                'CACACAGTGATC': 21
                }
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 994)
        self.assertEqual(barcode_counts, expected_barcode_counts)

    def testBarcodeSplitter(self):
        bc = bcs.BarcodeSplitter(in_file, out_dir, prefix + "barcodes.csv")
        num_read, num_written, barcode_counts = bc.run()
        expected_barcode_counts = {
                'CACACTTGAATC': 332,
                'CACCGATGTATC': 6,
                'CACATCACGATC': 7,
                'CACCAGATCATC': 281,
                'CACTTAGGCATC': 7,
                'CACGATCAGATC': 270,
                'CACGCCAATATC': 35,
                'CACTGACCAATC': 35,
                'CACACAGTGATC': 21
                }
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 994)
        self.assertEqual(barcode_counts, expected_barcode_counts)

    def testHardTrimmer(self):
        ht = htrim.HardTrimmer(in_file, out_dir + "ht.fastq", length=30)
        num_read, num_written = ht.run()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 1000)

    def testHardTrimmerParallel(self):
        ht = htrim.HardTrimmer(in_file, out_dir + "ht.fastq", length=30)
        num_read, num_written = ht.run_parallel()
        self.assertEqual(num_read, 1000)
        self.assertEqual(num_written, 1000)

    def testQualStats(self):
        qs = qstat.QualStats(in_file, qual_offset=33)
        self.assertTrue(qs.run())
#        num_read, num_written = qs.run()
#        self.assertEqual(num_read, 1000)
#        self.assertEqual(num_written, 1000)

    def testFastqToFasta(self):
        ftf = conv.FastqToFasta(in_file, out_dir + "fasta.fasta")
        self.assertTrue(ftf.run())
#        num_read, num_written = ftf.run()
#        self.assertEqual(num_read, 1000)
#        self.assertEqual(num_written, 1000)


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Tester))
    return suite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=3).run(suite())

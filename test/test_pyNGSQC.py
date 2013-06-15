#!/usr/bin/env python
# Copyright 2012 Kevin Murray
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pyngsqc as ngs
from pyngsqc import qualfilter as qfil
from pyngsqc import qualstats as qstat
from pyngsqc import qualtrimmer as qtrim
from pyngsqc import hardtrimmer as htrim
from pyngsqc import barcodesplitter as bcs
from pyngsqc import collapser as col
from pyngsqc import converter as conv
from test.data.expected import (
    EXPECTED_BARCODE_COUNTS,
    EXPECTED_QUALSTATS_POSITIONS,
)
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
        qf.run_parallel()
        self.assertEqual(qf.stats["reader"]["num_reads"], 1000)
        self.assertEqual(qf.stats["runner"]["num_reads"], 996)

    def testQualFilter(self):
        qf = qfil.QualFilter(
            in_file,
            out_dir + "qf.fastq",
            qual_threshold=20,
            qual_offset=33,
            pass_rate=0.9,
            max_Ns=-1,
        )
        qf.run()
        self.assertEqual(qf.stats["reader"]["num_reads"], 1000)
        self.assertEqual(qf.stats["writer"]["num_reads"], 996)

    def testCollapser(self):
        co = col.Collapser(
            in_file,
            out_dir + "col.fastq",
        )
        co.run()
        self.assertEqual(co.stats["reader"]["num_reads"], 1000)
        self.assertEqual(co.stats["writer"]["num_reads"], 992)

    def testQualTrimmerParallel(self):
        qt = qtrim.QualTrimmer(
            in_file,
            out_dir + "qt.fastq",
            qual_threshold=20,
            min_length=10,
            qual_offset=33
        )
        qt.run_parallel()
        self.assertEqual(qt.stats["reader"]["num_reads"], 1000)
        self.assertEqual(qt.stats["runner"]["num_reads"], 1000)

    def testQualTrimmer(self):
        qt = qtrim.QualTrimmer(
            in_file,
            out_dir + "qt.fastq",
            qual_threshold=20,
            min_length=10,
            qual_offset=33
        )
        qt.run()
        self.assertEqual(qt.stats["reader"]["num_reads"], 1000)
        self.assertEqual(qt.stats["writer"]["num_reads"], 1000)

    def testBarcodeSplitterParallel(self):
        bc = bcs.BarcodeSplitter(in_file, out_dir, prefix + "barcodes.csv")
        bc.run_parallel()
        self.assertEqual(bc.stats["runner"]["barcode_counts"],
                         EXPECTED_BARCODE_COUNTS)

    def testBarcodeSplitter(self):
        bc = bcs.BarcodeSplitter(in_file, out_dir, prefix + "barcodes.csv")
        bc.run()
        self.assertEqual(bc.stats["reader"]["num_reads"], 1000)
        self.assertEqual(bc.stats["writer"]["num_reads"], 994)
        self.assertEqual(bc.stats["writer"]["barcode_counts"],
                         EXPECTED_BARCODE_COUNTS)

    def testHardTrimmer(self):
        ht = htrim.HardTrimmer(in_file, out_dir + "ht.fastq", length=30)
        ht.run()
        self.assertEqual(ht.stats["reader"]["num_reads"], 1000)
        self.assertEqual(ht.stats["writer"]["num_reads"], 1000)

    def testHardTrimmerParallel(self):
        ht = htrim.HardTrimmer(in_file, out_dir + "ht.fastq", length=30)
        ht.run_parallel()
        self.assertEqual(ht.stats["reader"]["num_reads"], 1000)
        self.assertEqual(ht.stats["runner"]["num_reads"], 1000)

    def testQualStats(self):
        qs = qstat.QualStats(in_file, qual_offset=33)
        qs.run()
        self.assertEqual(qs.stats["reader"]["num_reads"], 1000)
        self.assertEqual(qs.stats["positions"], EXPECTED_QUALSTATS_POSITIONS)

    def testFastqToFasta(self):
        ftf = conv.FastqToFasta(in_file, out_dir + "fasta.fasta")
        ftf.run()
        self.assertEqual(ftf.stats["reader"]["num_reads"], 1000)
        self.assertEqual(ftf.stats["writer"]["num_reads"], 1000)

    def testConvertQualOffset(self):
        cpo = conv.ConvertQualOffset(
            in_file,
            out_dir + "phred64.fastq",
            in_qual_offset=33,
            out_qual_offset=64,
        )
        cpo.run()
        self.assertEqual(cpo.stats["reader"]["num_reads"], 1000)
        self.assertEqual(cpo.stats["writer"]["num_reads"], 1000)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=3)
    unittest.main(testRunner=runner)

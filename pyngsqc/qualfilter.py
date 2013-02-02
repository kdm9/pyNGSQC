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

import pyngsqc
from sys import stderr
import _parallel


class QualFilter(pyngsqc.QualBase):
    """
    Usage:
        QualityFilter(in_file_name, out_file_name, qual_threshold=15,
                  pass_rate=0.9, max_Ns=-1, qual_offset=64, append=False,
                  compression=pyngsqc.GUESS_COMPRESSION):
            in_file_name (str): path of input file, can be .fastq, .fastq.gz
                or .fastq.bz2
            out_file_name (str): path of output file, can be .fastq, .fastq.gz
                or .fastq.bz2
            qual_threshold (int): minimum "pass" phred score
            pass_rate (float): minimum fraction of bases which must be equal
                to or greater than qual_threshold
    """
    def __init__(
            self,
            # Inherited args
            in_file_name,
            out_file_name,
            # Local args
            # Inherited kwargs
            qual_offset=64,
            compression=pyngsqc.GUESS_COMPRESSION,
            deduplicate_header=True,
            verbose=False,
            print_summary=False,
            # Local kwargs
            qual_threshold=15,
            pass_rate=0.9,
            max_Ns=-1
            ):
        # Initialise base class
        super(QualFilter, self).__init__(
            in_file_name,
            out_file_name,
            qual_offset=qual_offset,
            compression=compression,
            deduplicate_header=deduplicate_header,
            verbose=verbose,
            print_summary=print_summary
            )
        # Initialise local variables
        self.pass_rate = float(pass_rate)
        self.qual_threshold = qual_threshold
        self.max_Ns = max_Ns

    def _passes_score_qc(self, read):
        qual = read[3]
        read_len = len(qual)
        low_scores = 0
        for p in qual:
            if pyngsqc.get_qual_from_phred(p, self.qual_offset) < self.qual_threshold:
                low_scores += 1
        this_pass_rate = 1.0 - float(low_scores) / float(read_len)
        if this_pass_rate <= self.pass_rate:
            return False
        else:
            return True

    def _passes_qc(self, read):
        if self.max_Ns != -1 and \
         int(pyngsqc.num_Ns_in_read(read)) > self.max_Ns:
            self.num_bad_reads += 1
            return False
        elif not self._passes_score_qc(read):
            self.num_bad_reads += 1
            return False
        else:
            self.num_good_reads += 1
            return True

    def _print_summary(self):
        stderr.write("QC check finished:\n")
        stderr.write("Processed %i reads\n" % self.num_reads)
        stderr.write(
            "\t%i sequences passed QC, wrote them to %s\n" %
            (self.num_good_reads, self.out_file_name)
            )
        stderr.write(
            "\t%i sequences failed QC, and were ignored\n" %
            self.num_bad_reads,
            )

    def filter_read(self, read):
        if self._passes_qc(read):
            return read
        else:
            return []

    def run(self):
        for read in self.reader:
            if self.filter_read(read):
                self.writer.write(read)
        self.num_reads += self.reader.num_reads
        if self.print_summary:
            self._print_summary()
        self.reader.close()
        return True

    def run_parallel(self):
        task_args = (
            self.pass_rate,
            self.qual_threshold,
            self.qual_offset,
            self.max_Ns
            )

        runner = _parallel.ParallelRunner(
            QualFilterTask,
            self.reader,
            self.writer,
            task_args
            )
        runner.run()
        self.num_good_reads = runner.writer.num_reads
        self.num_reads = runner.num_reads
        if self.print_summary:
            self._print_summary()
        return True


class QualFilterTask(QualFilter):

    def __init__(
                 self,
                 read,
                 qual_threshold,
                 pass_rate,
                 max_Ns,
                 qual_offset,
                ):
        self.read = read
        self.pass_rate = float(pass_rate)
        self.qual_threshold = qual_threshold
        self.qual_offset = qual_offset
        self.max_Ns = max_Ns

    def __call__(self):
        return self.filter_read(self.read)

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

from sys import stderr
import pyngsqc
import _paralell


class QualTrimmer(pyngsqc.QualBase):

    def __init__(
            self,
            # Inherited args
            in_file_name,
            out_file_name,
            # Local args
            # Inherited kwargs
            qual_offset=pyngsqc.DEFAULT_QUAL_OFFSET,
            qual_threshold=pyngsqc.DEFAULT_QUAL_THRESHOLD,
            compression=pyngsqc.GUESS_COMPRESSION,
            deduplicate_header=True,
            verbose=False,
            # Local kwargs
            min_length=15,
            remove_trailing_Ns=False,
            ):
        # Initialise base class
        super(QualTrimmer, self).__init__(
            in_file_name,
            out_file_name,
            qual_offset=qual_offset,
            qual_threshold=qual_threshold,
            compression=compression,
            deduplicate_header=deduplicate_header,
            verbose=verbose
            )
        # Initialise local variables
        self.remove_trailing_Ns = remove_trailing_Ns
        self.min_length = min_length

    def trim_read(self, read):
        for iii in reversed(xrange(len(read[1]))):
            if self.remove_trailing_Ns and read[1][iii].upper() == "N":
                # Remove last base and score
                read[1] = read[1][:-1]  # Base
                read[3] = read[3][:-1]  # Phred Score
            elif self._get_qual_from_phred(read[3][iii]) < \
             self.qual_threshold:
                read[1] = read[1][:-1]  # Base
                read[3] = read[3][:-1]  # Phred Score
        if len(read[1]) < self.min_length:
            return None
        else:
            return read

    def print_summary(self):
        stderr.write("QualTrimmer finished:\n")
        stderr.write(
                      "\t%i sequences passed QC, wrote them to %s\n" %
                      (self.num_good_reads, self.out_file_name)
                    )
        stderr.write(
                      "\t%i sequences failed QC, and were ignored\n" %
                      self.num_bad_reads,
                    )

    def run(self):
        for read in self.reader:
            self.num_reads += 1
            read = self.trim_read(read)
            if read is not None:
                self.num_good_reads += 1
                self.writer.write(read)
            else:
                self.num_bad_reads += 1
        self.print_summary()
        return True

    def run_paralell(self):
        runner = _paralell.ParalellRunner(
            QualTrimmerTask,
            self.reader,
            self.writer,
            (
                self.qual_threshold,
                self.qual_offset,
                self.min_length,
                self.remove_trailing_Ns,
                )
            )
        runner.run()
        self.num_good_reads = runner.writer.num_reads
        self.num_reads = runner.num_reads
        self.print_summary()
        return True


class QualTrimmerTask(QualTrimmer):

    def __init__(
                 self,
                 read,
                 qual_threshold,
                 qual_offset,
                 min_length,
                 remove_trailing_Ns,
                ):
        self.read = read
        self.qual_threshold = qual_threshold
        self.qual_offset = qual_offset
        self.min_length = min_length
        self.remove_trailing_Ns = remove_trailing_Ns

    def __call__(self):
        return self.trim_read(self.read)

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
import paralellNGS
import pyNGSQC


class HardTrimmer(pyNGSQC.NGSQC):

    def __init__(
                 self,
                 in_file_name,
                 out_file_name,
                 length=15,
                 verbose=False,
                 compression=pyNGSQC.GUESS_COMPRESSION
                ):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.reader = pyNGSQC.FastqReader(
                                           self.in_file_name,
                                           compression=compression
                                         )
        self.writer = pyNGSQC.FastqWriter(
                                           self.out_file_name,
                                           compression=compression
                                          )
        self.verbose = verbose
        self.length = length
        self.num_reads = 0

    def trim_read(self, read):
        read[1] = read[1][:self.length]  # Base
        read[3] = read[3][:self.length]  # Phred Score
        return read

    def print_summary(self):
        stderr.write("HardTrimmer finished:\n")
        stderr.write(
                      "\t%i sequences were trimmed to %i wrote them to %s\n" %
                      (self.num_reads, self.length, self.out_file_name)
                    )

    def run(self):
        for read in self.reader:
            read = self.trim_read(read)
            if read is not None:
                self.num_reads += 1
                self.writer.write(read)
        self.print_summary()
        return True

    def run_paralell(self):
        runner = paralellNGS.ParalellRunner(
            HardTrimmerTask,
            self.reader,
            self.writer,
            (self.length,)
            )
        runner.run()
        self.num_reads = runner.num_reads
        self.print_summary()
        return True


class HardTrimmerTask(HardTrimmer):

    def __init__(self, read, length):
        self.read = read
        self.length = length

    def __call__(self):
        return self.trim_read(self.read)

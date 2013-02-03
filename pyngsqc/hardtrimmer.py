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
import _parallel
import pyngsqc


class HardTrimmer(pyngsqc.Base):
    """
    """
    def __init__(
            self,
            in_file_name,
            out_file_name,
            length=15,
            verbose=False,
            compression=pyngsqc.GUESS_COMPRESSION,
            print_summary=False
            ):
        super(HardTrimmer, self).__init__(
            in_file_name,
            out_file_name,
            verbose=verbose,
            compression=compression,
            print_summary=print_summary
            )
        self.length = length

    def trim_read(self, read):
        read[1] = read[1][:self.length]  # Base
        read[3] = read[3][:self.length]  # Phred Score
        return read

    def _print_summary(self):
        stderr.write("HardTrimmer finished:\n")
        stderr.write(
            "\t%i sequences were trimmed to %i wrote them to %s\n" %
            (self.writer.stats["num_reads"], self.length, self.out_file_name)
            )

    def run(self):
        for read in self.reader:
            read = self.trim_read(read)
            if len(read) == 4:
                self.writer.write(read)
        if self.print_summary:
            self._print_summary()
        return (
                self.reader.stats["num_reads"],
                self.writer.stats["num_reads"]
                )

    def run_parallel(self):
        runner = _parallel.ParallelRunner(
            HardTrimmerTask,
            self.reader,
            self.writer,
            (self.length,)
            )
        runner.run()
        if self.print_summary:
            self._print_summary()
        return (
                self.reader.stats["num_reads"],
                runner.stats["num_reads"]
                )


class HardTrimmerTask(HardTrimmer):

    def __init__(self, read, length):
        self.read = read
        self.length = length

    def __call__(self):
        return self.trim_read(self.read)

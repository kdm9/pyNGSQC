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

import pyngsqc
from pyngsqc import _parallel
import sys
from tempfile import NamedTemporaryFile as namedtmp
import os


class Collapser(pyngsqc.Base):
    # MAX_FILE_SIZE = 2 << 30  # 2GB

    def __init__(
            self,
            # Inherited args
            in_file_name,
            out_file_name,
            # Local args
            # Inherited kwargs
            compression="guess",
            deduplicate_header=True,
            verbose=False,
            print_summary=False,
            # Local kwargs
            key_length=5,
            tmp_dir=None
    ):
        # Initialise base class
        super(Collapser, self).__init__(
            in_file_name,
            out_file_name,
            compression=compression,
            deduplicate_header=deduplicate_header,
            verbose=verbose,
            print_summary=print_summary
        )
        # Initialise local variables
        self.tmp_dir = tmp_dir
        self.key_length = key_length
        self.keys = []
        self.tmp_file_names = {}
        self.file_sizes = {}

    def _split_files(self):
        """
        Splits input file into subfiles based on their "key", or first
         key_length bases, to allow for a memory efficient sort
        """
        for read in self.reader:
            key = read[1][:self.key_length]
            if key in self.tmp_file_names:
                fh = open(self.tmp_file_names[key], "a")
            else:
                self.keys.append(key)
                # If in keys, file handle should exist
                fh = namedtmp(
                    mode="w",
                    dir=self.tmp_dir,
                    prefix=key + "_",
                    delete=False
                )
                file_name = fh.name
                self.tmp_file_names[key] = file_name
            read_str = "\n".join(read) + "\n"
            fh.write(read_str)
            fh.close()

        # get file size
        for key in self.tmp_file_names:
            fh = open(self.tmp_file_names[key], "rb")
            fh.seek(0, 2)  # go the end of the file
            this_file_size = fh.tell()  # and get its size
            self.file_sizes[key] = this_file_size
            fh.close

    def _read_to_tuple(self, read):
        return (read[1], read[0], read[3], read[2])

    def _tuple_to_read(self, read_tuple):
        read = []
        read.append(read_tuple[1])
        read.append(read_tuple[0])
        read.append(read_tuple[3])
        read.append(read_tuple[2])
        return read

    def _colapse(self):
        for key in sorted(self.keys):
            these_reads = []
            file_name = self.tmp_file_names[key]
            reader = pyngsqc.FastqReader(file_name)
            for read in reader:
                these_reads.append(self._read_to_tuple(read))
            these_reads.sort()
            reader.close()
            last_read_seq = ""
            for read_tuple in these_reads:
                if read_tuple[0] != last_read_seq:
                    last_read_seq = read_tuple[0]
                    self.writer.write(self._tuple_to_read(read_tuple))
        for file_name in list(self.tmp_file_names.values()):
            os.remove(file_name)

    def _print_summary(self):
        sys.stderr.write("Collapser finished\n")
        sys.stderr.write("\tAnalysed %i reads\n" %
                         self.stats["reader"]["num_reads"])
        sys.stderr.write("\tFound %i unique reads\n" %
                         self.stats["writer"]["num_reads"])
        sys.stderr.write(
            "\tRemoved %i non-unique reads\n" %
            self.stats["reader"]["num_reads"] - self.stats[
            "writer"]["num_reads"]
        )

    def run(self):
        self._split_files()
        self._colapse()

        self.stats["reader"] = self.reader.stats
        self.stats["writer"] = self.writer.stats
        if self.print_summary:
            self._print_summary()

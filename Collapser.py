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

import pyNGSQC
import paralellNGS
import sys
from tempfile import NamedTemporaryFile as namedtmp
import os


class Collapser(pyNGSQC.Base):
    #MAX_FILE_SIZE = 2 << 30  # 2GB

    def __init__(
            self,
            in_file_name,
            out_file_name,
            key_length=5,
            tmp_dir=None,
            compression=pyNGSQC.GUESS_COMPRESSION,
            verbose=False
            ):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.reader = pyNGSQC.FastqReader(
            self.in_file_name,
            compression=compression
            )
        self.tmp_dir = tmp_dir
        self.verbose = verbose
        self.key_length = key_length
        self.keys = []
        self.num_reads = 0L
        self.num_non_unique_reads = 0
        self.num_unique_reads = 0
        self.tmp_file_names = {}
        self.file_sizes = {}

    def _split_files(self):

        for read in self.reader:
            key = read[1][:self.key_length]
            # None means guess tmp dir
            if key in self.tmp_file_names:
                fh = open(self.tmp_file_names[key], "ab")
            else:
                self.keys.append(key)
                # If in keys, file handle should exist
                fh = namedtmp(
                    mode="wb",
                    dir=self.tmp_dir,
                    prefix=key + "_",
                    delete=False
                    )
                file_name = fh.name
                self.tmp_file_names[key] = file_name
            read_str = "\n".join(read)
            fh.write(read_str + "\n")
            fh.close()

        # get file size
        for key in self.tmp_file_names:
            fh = open(self.tmp_file_names[key], "rb")
            fh.seek(0, 2)
            this_file_size = fh.tell()
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
        sorted_writer = pyNGSQC.FastqWriter(self.out_file_name)
        for key in sorted(self.keys):
            these_reads = []
            file_name = self.tmp_file_names[key]

            reader = pyNGSQC.FastqReader(file_name)
            for read in reader:
                these_reads.append(self._read_to_tuple(read))
            these_reads.sort()
            self.num_reads += reader.num_reads
            reader.close()

            last_read_seq = ""
            for read_tuple in these_reads:

                if read_tuple[0] != last_read_seq:
                    last_read_seq = read_tuple[0]
                    self.num_unique_reads += 1
                    sorted_writer.write(self._tuple_to_read(read_tuple))
                else:
                    self.num_non_unique_reads += 1
        for file_name in self.tmp_file_names.values():
            os.remove(file_name)

    def print_summary(self):
        sys.stderr.write("Collapser finished\n")
        sys.stderr.write("\tAnalysed %i reads\n" % self.num_reads)
        sys.stderr.write("\tFound %i unique reads\n" % self.num_unique_reads)
        sys.stderr.write("\tRemoved %i non-unique reads\n" %
         self.num_non_unique_reads)

    def run(self):
        self._split_files()
        print "Files Split, sizes:"
        for key, size in self.file_sizes:
            print "\t%s %r" % (key, size)
        self._colapse()
        print "collapsed"
        self.print_summary()
        return True

    def run_paralell(self):
        # We don't bother paralellising spliting of files, as it is mostly IO,
        # so would be faily pointless
        self._split_files()

        sorted_writer = pyNGSQC.FastqWriter(self.out_file_name)
        files = pyNGSQC.dict_to_tuples(self.tmp_file_names)
        runner = paralellNGS.ParalellRunner(
            CollapserTask,
            files,
            sorted_writer
            )
        runner.run()
        self.barcode_counts = runner.writer.barcode_counts
        self.num_reads = runner.num_reads
        self.print_summary()
        return True
        self._sort()
        self.print_summary()
        return True


class CollapserTask(Collapser):

    def __init__(self, file_tuple):
        self.file_tuple = file_tuple
        self.num_reads = 0L
        self.num_non_unique_reads = 0
        self.num_unique_reads = 0

    def call(self):
        key, file_name = self.file_tuple
        these_reads = []
        these_unique_reads = []

        reader = pyNGSQC.FastqReader(file_name)
        for read in reader:
            these_reads.append(self._read_to_tuple(read))
        these_reads.sort()
        self.num_reads += reader.num_reads
        reader.close()

        last_read_seq = ""
        for read_tuple in these_reads:

            if read_tuple[0] != last_read_seq:
                last_read_seq = read_tuple[0]
                self.num_unique_reads += 1
                for line in self._tuple_to_read(read_tuple):
                    these_unique_reads.append(line)
            else:
                self.num_non_unique_reads += 1
        os.remove(file_name)

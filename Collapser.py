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


class Collapser(pyNGSQC.NGSQC):
    MAX_FILE_SIZE = 2 << 30  # 2GB

    def __init__(
                 self,
                 in_file_name,
                 out_file_name,
                 key_length=4,
                 compression=pyNGSQC.GUESS_COMPRESSION,
                 verbose=False
                ):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.reader = pyNGSQC.FastqReader(
                                           self.in_file_name,
                                           compression=compression
                                         )
        self.verbose = verbose
        self.key_length = key_length
        self.keys = []
        self.num_reads = 0L
        self.tmp_file_name = "%s.tmp"

    def _split_files(self):
        file_sizes = {}
        sub_file_handles = {}
        for read in self.reader:
            key = read[1][:self.key_length]
            file_name = self.tmp_file_name % key
            if key not in self.keys:
                self.keys.append(key)
                # If in keys, file handle should exist
                sub_file_handles[key] = open(file_name, "wb")
            read_str = "\n".join(read)
            sub_file_handles[key].write(read_str + "\n")

        # get file size
        for key in self.keys:
            sub_file_handles[key].seek(0, 2)
            this_file_size = sub_file_handles[key].tell()
            file_sizes[key] = this_file_size
            sub_file_handles[key].close

#        # Not sure if this works, needs debugging
#        while max(file_sizes.values()) > self.MAX_FILE_SIZE:
#            self.key_length += 1
#            for parent_key in self.keys:
#                this_read = []
#                parent_reader = pyNGSQC.FastqReader(in_file_name)
#                sub_file_handles[parent_key].close()
#                for line in self.sub_files[parent_key]:
#                    this_read.append(line)
#                    if len(this_read) == 4:
#                        this_key = this_read[1][:self.key_length]
#                        file_name = "collaper_sort_%s" % this_key
#                        if this_key not in self.sub_files:
#                            self.sub_files[this_key] = open(file_name, "wb")
#                        read_str = "\n".join(this_read)
#                        self.sub_files[this_key].write(read_str + "\n")
#                sub_file_handles[key].seek(0)
#                this_file_size = sub_file_handles[key].tell()
#                file_sizes[key] = this_file_size
#                del file_sizes[parent_key]
#                sub_file_handles[parent_key].close()
#                del sub_file_handles[parent_key]

    def _read_to_tuple(self, read):
        return (read[1], read[0], read[3], read[2])

    def _tuple_to_read(self, read_tuple):
        read = []
        read.append(read_tuple[1])
        read.append(read_tuple[0])
        read.append(read_tuple[3])
        read.append(read_tuple[2])
        return read

    def _sort(self):
        for key in sorted(self.keys):
            these_reads = []
            file_name = self.tmp_file_name % key

            reader = pyNGSQC.FastqReader(file_name)
            for read in reader:
                these_reads.append(self._read_to_tuple(read))
            these_reads.sort()
            reader.close()

            sorted_writer = pyNGSQC.FastqWriter("sorted_" + file_name)
            last_read_seq = ""
            for read_tuple in these_reads:
                if read_tuple[0] != last_read_seq:
                    sorted_writer.write(self._tuple_to_read(read_tuple))
                    last_read_seq = read_tuple[0]

    def print_summary(self):
        pass

    def run(self):
        self._split_files()
        self._sort()
        return True
fldr = "/home/kevin/workspace/"  # /var/ws/borevitz/
cl = Collapser(
               fldr + "in.fastq",  # pel_1.small.fastq
               fldr + "test.out.fastq"
              )
cl.run()

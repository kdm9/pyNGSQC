#from sys import stderr
import pyNGSQC
#import os
#import tempfile


class Collapser(pyNGSQC.NGSQC):
    MAX_FILE_SIZE = 2 << (3 * 10)  # 2GB

    def __init__(
                 self,
                 in_file_name,
                 out_file_name,
                 key_length=2,
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
        self.sub_files = {}
        self.num_reads = 0L

    def split_files(self):
        file_sizes = {}

        for read in self.reader:
            this_key = read[1][:self.key_length]
            file_name = "collaper_sort_%s" % this_key
            if this_key not in self.sub_files:
                self.sub_files[this_key] = open(file_name, "wb")
            read_str = "\n".join(read)
            self.sub_files[this_key].write(read_str + "\n")

        for key in self.sub_files:
            # get file size
            self.sub_files[key].seek(0, 2)
            this_file_size = self.sub_files[key].tell()
            file_sizes[key] = this_file_size

        # Not sure if this works, needs debugging
        while max(file_sizes.values()) > self.MAX_FILE_SIZE:
            self.key_length += 1
            for parent_key in self.sub_files:
                this_read = []
                self.sub_files[parent_key].seek(0)
                for line in self.sub_files[parent_key]:
                    this_read.append(line)
                    if len(this_read) == 4:
                        this_key = this_read[1][:self.key_length]
                        file_name = "collaper_sort_%s" % this_key
                        if this_key not in self.sub_files:
                            self.sub_files[this_key] = open(file_name, "wb")
                        read_str = "\n".join(this_read)
                        self.sub_files[this_key].write(read_str + "\n")
                self.sub_files[key].seek(0)
                this_file_size = self.sub_files[key].tell()
                file_sizes[key] = this_file_size
                del file_sizes[parent_key]
                self.sub_files[parent_key].close()
                del self.sub_files[parent_key]

    def sort(self):
        sub_file_readers = {}
        sorted_file_writers = {}
        for key in self.sub_files:
            self.sub_files[key].close()
            file_name = "collaper_sort_%s" % key
            sub_file_readers[key] = pyNGSQC.FastqReader(file_name)
            #del self.sub_files[key]

        for key in sub_file_readers:
            print key
            out_file_name = "collaper_sorted_%s" % key
            in_file_name = "collaper_sort_%s" % key
            these_reads = []
            sorted_file_writers[key] = pyNGSQC.FastqWriter(out_file_name)
            this_reader = pyNGSQC.FastqReader(in_file_name)
            for read in this_reader:
                these_reads.append((read[1], read[0]))
            these_reads.sort()
            for pair in these_reads:
                print pair

    def print_summary(self):
        pass

    def run(self):
        return True

cl = Collapser(
               "/var/ws/borevitz/pel_1.small.fastq",
               "/var/ws/borevitz/test.out.fastq"
              )
cl.split_files()
cl.sort()

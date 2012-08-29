#from sys import stderr
import pyNGSQC
#from os import path
import tempfile


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
            if this_key not in self.sub_files:
                self.sub_files[this_key] = tempfile.TemporaryFile()
            read_str = "\n".join(read)
            self.sub_files[this_key].write(read_str + "\n")

        for key in self.sub_files:
            # get file size
            self.sub_files[key].seek(0,2)
            this_file_size = self.sub_files[key].tell()
            file_sizes[key] = this_file_size

        while max(file_sizes.values()) > self.MAX_FILE_SIZE:
            self.key_length += 1
            for parent_key in self.sub_files:
                this_read = []
                self.sub_files[parent_key].seek(0)
                for line in self.sub_files[parent_key]:
                    this_read.append(line)
                    if len(this_read) == 4:
                        this_key = this_read[1][:self.key_length]
                        if this_key not in self.sub_files:
                            self.sub_files[this_key] = tempfile.TemporaryFile()
                        read_str = "\n".join(this_read)
                        self.sub_files[this_key].write(read_str + "\n")
                self.sub_files[key].seek(0)
                this_file_size = self.sub_files[key].tell()
                file_sizes[key] = this_file_size
                del file_sizes[parent_key]
                self.sub_files[parent_key].close()
                del self.sub_files[parent_key]
        print file_sizes
    def print_summary(self):
        pass

    def run(self):
        return True

cl = Collapser(
               "/var/ws/borevitz/test.fastq",
               "/var/ws/borevitz/test.out.fastq"
              )
cl.split_files()

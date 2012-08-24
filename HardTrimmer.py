from sys import stderr
import pyNGSQC


class HardTrimmer(pyNGSQC.NGSQC):
    """
    """
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

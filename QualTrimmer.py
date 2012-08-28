from sys import stderr
import pyNGSQC


class QualTrimmer(pyNGSQC.NGSQC):
    def __init__(
                 self,
                 in_file_name,
                 out_file_name,
                 qual_threshold=20,
                 qual_offset=64,
                 min_length=15,
                 remove_trailing_Ns=False,
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
        self.qual_threshold = qual_threshold
        self.remove_trailing_Ns = remove_trailing_Ns
        self.qual_offset = qual_offset
        self.min_length = min_length
        self.num_reads = 0
        self.num_good_reads = 0
        self.num_bad_reads = 0

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

import pyNGSQC
from sys import stderr


class QualStats(pyNGSQC.NGSQC):
    """
    Usage:
    """

    def __init__(self, in_file_name, out_file_name, qual_threshold=15,
                  pass_rate=0.9, max_Ns=-1, qual_offset=64,
                  compression=pyNGSQC.GUESS_COMPRESSION):
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
        self.pass_rate = float(pass_rate)
        self.qual_threshold = qual_threshold
        self.qual_offset = qual_offset
        self.max_Ns = max_Ns
        self.num_reads = 0
        self.num_good_reads = 0
        self.num_bad_reads = 0

    def _passes_score_qc(self, read):
        qual = read[3]
        read_len = len(qual)
        low_scores = 0
        for p in qual:
            if self._get_qual_from_phred(p) < self.qual_threshold:
                low_scores += 1
        this_pass_rate = 1.0 - float(low_scores) / float(read_len)
        if this_pass_rate <= self.pass_rate:
            return False
        else:
            return True

    def filter_read(self, read):
        if self.max_Ns != -1 and\
         int(self._num_Ns_in_read(read)) > self.max_Ns:
            self.num_bad_reads += 1
            return False
        elif not self._passes_score_qc(read):
            self.num_bad_reads += 1
            return False
        else:
            self.num_good_reads += 1
            return True

    def print_summary(self):
        stderr.write("QC check finished:\n")
        stderr.write("Processed %i reads\n" % self.num_reads)
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
            if self.filter_read(read):
                self.writer.write(read)
        self.print_summary()
        return True

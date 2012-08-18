from sys import stderr 
import gzip
import bz2
import pyNGSQC

class QualityTrimmer(pyNGSQC.NGSQC):
    def __init__(
                 self,
                 in_file_name,
                 out_file_name,
                 append=False,
                 threshold=20,
                 qual_offset=64,
                 min_length=15,
                 remove_trailing_Ns=True,
                 compression=pyNGSQC.GUESS_COMPRESSION
                ):
        
        super(QualityFilter, self).__init__(
                                            in_file_name,
                                            out_file_name,
                                            qual_offset=qual_offset,
                                            append=append,
                                            compression=compression
                                           )
        self.remove_trailing_Ns = remove_trailing_Ns
        self.min_quality = qual_threshold
        self.num_good_reads = 0
        self.num_bad_reads = 0
        self.num_records = -1
    
    def trim_read(self, read):
        for iii in reversed(xrange(len(read[1]))):
            if read[1][i].upper() == "N" and self.remove_trailing_Ns:
                # Remove last base and score
                read[1].pop() # Base
                read[3].pop() # Phred Score
            if self.get_qual_from_phred(read[3][i]) < self.min_quality:
                read[1].pop()
                read[3].pop()
        if len(read[1]) < self.min_length:
            return None
        else:
            return read
    
    def print_summary(self):
        stderr.write("QC check finished:\n")
        stderr.write(
                      "\t%i sequences passed QC, wrote them to %s\n" % 
                      (self.num_good_reads,  self.out_file_name)
                    )
        stderr.write(
                      "\t%i sequences failed QC, and were ignored\n" %
                      self.num_bad_reads,
                    )
    
    def process_read(self, this_read):
        if self.trim_read(this_read) is not None:
            self._write_good_read(this_read)

if __name__ == "__main__":
    import sys
    qf = QualityTrimmer(sys.argv[1],
                       sys.argv[2],
                       compression=int(sys.argv[3]),
                       qual_offset=33)
    qf.run(qf.process_read, qf.print_summary)


from sys import stderr 
import gzip
import bz2
import pyNGSQC

class QualityFilter(pyNGSQC.NGSQC):
    """
    Usage:
        QualityFilter(in_file_name, out_file_name, qual_threshold=15,
                  pass_rate=0.9, max_Ns=-1, qual_offset=64, append=False,
                  compression=pyNGSQC.GUESS_COMPRESSION):
            in_file_name (str): path of input file, can be .fastq, .fastq.gz
                or .fastq.bz2
            out_file_name (str): path of output file, can be .fastq, .fastq.gz
                or .fastq.bz2
            qual_threshold (int): minimum "pass" phred score
            pass_rate (float): minimum fraction of 
    """
    def __init__(self, in_file_name, out_file_name, qual_threshold=15,
                  pass_rate=0.9, max_Ns=-1, qual_offset=64, append=False,
                  compression=pyNGSQC.GUESS_COMPRESSION):
        
        super(QualityFilter, self).__init__(in_file_name, out_file_name,
                  qual_offset=qual_offset, append=append,
                  compression=compression)
        
        self.pass_rate = float(pass_rate)
        self.qual_threshold = qual_threshold
        self.max_Ns = max_Ns
    
    def _has_Ns(self, read):
        #seq_header, seq, qual_header, qual = read
        seq = read[1]
        seq = seq.upper()
        if "N" not in seq:
            return False
        else:
            return seq.count('N')

    def _has_low_score(self,read):
        qual = read[3]
        read_len = len(qual)
        low_scores = 0
        for p in qual:
            if self.get_qual_from_phred(p) < self.qual_threshold:
                low_scores += 1
        this_pass_rate = 1.0 - float(low_scores)/float(read_len)
        if this_pass_rate < self.pass_rate:
            return False
        else:
            return True
            
    def filter_read(self, read):
        # the and max_Ns >= 0 skips the N-check if not required
        if self.max_Ns >= 0 and int(self._has_Ns(read)) > self.max_Ns :
            self.num_bad_reads += 1
            return False
        elif self._has_low_score(read):
            self.num_bad_reads += 1
            return False
        else:
            self.num_good_reads += 1
            return True
    
    def process_read(self, this_read):
        if self.filter_read(this_read):
                    self._write_good_read(this_read)
    
    def run(self):
        self._run(self.process_read, self.print_summary)

if __name__ == "__main__":
    import sys
    qf = QualityFilter(sys.argv[1],
                       sys.argv[2],
                       compression=int(sys.argv[3]),
                       pass_rate=1.0, 
                       qual_offset=33)
    qf.run()

from sys import stderr 
NO_ZIP = 0
GZIPED = 1
BZIP2ED = 2
class QualityFilter(object):
    def __init__(self, in_file_name, out_file_name, min_quality=15,
                  max_Ns=0, qual_offset=64, append=False, compression=NO_ZIP):
        self.compression = compression
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.min_quality = min_quality
        self.max_Ns = max_Ns
        self.qual_offset = qual_offset
        self.num_good_reads = 0
        self.num_bad_reads = 0
        self.num_records = -1
        if self.compression == NO_ZIP:
            if append:
                write_mode = "ab"
            else:
                write_mode = "wb"
            self.in_file = open(self.in_file_name, "rb")
            self.out_file = open(self.out_file_name, write_mode)
        elif self.compression == GZIPED:
            self.in_file = ""
            self.out_file = ""
        elif self.compression == BZIP2ED:
            self.in_file = ""
            self.out_file = ""
        else:
            raise ValueError(self.compression)
    
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
        for p in qual:
            if self.get_qual_from_phred(p) < self.min_quality:
                return False
        return True
            
    def qc_read(self, read):
        #print read[1], int(self._has_Ns(read))
        if int(self._has_Ns(read)) > self.max_Ns:
            self.num_bad_reads += 1
            return False
#        elif self._has_low_score(read):
#            self.num_bad_reads += 1
#            return False
        else:
            self.num_good_reads += 1
            return True
    
    def get_qual_from_phred(self, phred):
        qual = ord(phred) - self.qual_offset
        if qual <0:
            raise ValueError(
                              "Invalid quality score %i from phred %s" %
                              (qual, phred)
                            )
        else:
            return qual
    
    
    def _write_good_read(self, read):
        for line in read:
            self.out_file.write(line + "\n")
    
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
    
    def run(self):
        this_read = []
        
        for line in self.in_file:
            if self.num_records < 0:
                if line[0] != "@":
                    continue
                else:
                    self.num_records += 1
            this_read.append(line)
            if len(this_read) == 4:
                if self.qc_read(this_read):
                    self._write_good_read(this_read)
                this_read = []
                self.num_records += 1
        self.print_summary()
        if self.compression == NO_ZIP:
            self.in_file.close()
            self.out_file.close()

if __name__ == "__main__":
    import sys
    qf = QualityFilter(sys.argv[1], sys.argv[2])
    qf.run()

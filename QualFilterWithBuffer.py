from sys import stderr 
MAX_NS = 0
class QualityFilter(object):
    def __init__(self, in_file_name, out_file_name, min_quality=15,
                  qual_offset=64, buffer_size=10000000, append=False):
        # buffer_size = 1 is approx 100b, so 10000000 lines = 1Mb
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.min_quality = min_quality
        self.qual_offset = qual_offset
        self.buffer_size = buffer_size
        self.buffer = []
        self.append = append
        self.num_good_reads = 0
        self.num_bad_reads = 0
        if self.append:
            write_mode = "ab"
        else:
            write_mode = "wb"
        self.in_file = open(self.in_file_name, "rb")
        self.out_file = open(self.out_file_name, write_mode)

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
        if int(self._has_Ns(read)) > MAX_NS:
            self.num_bad_reads += 1
            return False
        elif self._has_low_score(read):
            self.num_bad_reads += 1
            return False
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
    
    def _load_buffer(self):
        while len(self.buffer) < self.buffer_size:
            line = self.in_file.readline()
            if line == "":
                return False
            else:
                self.buffer.append(line)
        long_enough = len(self.buffer) >= 4
        if not long_enough:
            print self.buffer
        return long_enough
    
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
        #loop thru file, using buffer
        while self._load_buffer():
            #if there is at least 
            while len(self.buffer) >= 4:
                this_read = []
                #first_line = 
                # If the first letter of the first line in the buffer is not @,
                # delete the first line of the buffer. this should only happen
                # at the start of the file as if we hit the end of the buffer
                # half-way through
                while self.buffer[0][0] != "@":
                    self.buffer.pop(0)
                    if len(self.buffer) < 4:
                        self.print_summary()
                        return
                #from here, we should be able to navigate by 4's 
                for r in xrange(4):
                    this_read.append(self.buffer.pop(0).rstrip())
                if self.qc_read(this_read):
                    self._write_good_read(this_read)
        self.print_summary()
        self.in_file.close()
        self.out_file.close()

if __name__ == "__main__":
    qf = QualityFilter("./in.fastq", "./out.fq")
    qf.run()

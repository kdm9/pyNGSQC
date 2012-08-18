import gzip
import bz2
import pyNGSQC
import os.path

GUESS_COMPRESSION = 0
NO_COMPRESSION = 1
GZIPPED = 2
BZIPP2ED = 3

BARCODE_FORWARD_ONLY = 0

SNIFF_CSV_DIALECT = 0
class NGSQC(object):
    
    def __init__(self, in_file_name, out_file_name, qual_offset=64,
                  append=False, compression=GUESS_COMPRESSION):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.qual_offset = qual_offset
        self.num_good_reads = 0
        self.num_bad_reads = 0
        self.num_records = -1
        
    # Sub-methods:
    
    def _num_Ns_in_read(self, read):
        #seq_header, seq, qual_header, qual = read
        seq = read[1]
        seq = seq.upper()
        if "N" not in seq:
            return False
        else:
            return seq.count('N')

    
    def _get_qual_from_phred(self, phred):
        qual = ord(phred) - self.qual_offset
        if qual < 0:
            raise ValueError(
                              "Invalid quality score %i from phred %s (ord %i)"\
                              % (qual, phred, ord(phred))
                            )
        return qual
    
    
    # General Methods:
    
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
    
    
    # Inherited methods, used only within children
    def _run(self, read_method, final_method):
        this_read = []
        for line in self.in_file:
            line = line.strip()
            if self.num_records < 0:
                if line[0] != "@":
                    continue
                else:
                    self.num_records += 1
            this_read.append(line)
            if len(this_read) == 4:
                read_method(this_read)
                this_read = []
                self.num_records += 1
        final_method()
        self.in_file.close()
        self.out_file.close()

class FastqIO(object):
    READ = 0
    WRITE = 1
    #APPEND = 2 # APPEND NOT SUPPORTED YET
    def __init__(self, file_name, mode=READ, compression=GUESS_COMPRESSION):
        self.file_name = file_name
        if mode == self.READ or mode == self.WRITE:
            self.mode = mode
        else:
            raise ValueError("%i is not a valid IO mode" % mode)
        if compression == GUESS_COMPRESSION:
            self.compression = self._guess_compression(self.file_name)
        elif compression == NO_COMPRESSION or compression == GZIPPED or \
         compression == BZIPP2ED:
            self.compression = compression
        else:
            raise ValueError("%i is not a valid compression mode" % compression)
    
    def _guess_compression(self, file_name):
        path, ext = os.path.splitext(file_name)
        gz_exts = [".gz", ".gzip"]
        bz2_exts = [".bz", ".bz2", ".bzip2"]
        if ext in gz_exts:
            return GZIPPED
        elif ext in bz2_exts:
            return BZIPP2ED
        else:
            return NO_COMPRESSION
    
    def get(self):
        if self.compression == NO_COMPRESSION:
            return _get_plaintext()
        elif self.compression == GZIPPED:
            return _get_gzip()
        elif self.compression == BZIPP2ED:
            return _get_bzip2()
        else:
            raise ValueError(self.compression)
    
    def _get_plaintext(self):
        if mode == self.READ:
            self.in_file = open(self.file_name, "rb")
        elif mode == self.WRITE:
            self.out_file = open(self.file_name, "wb")
        
    
    def _get_gzip(self):
        if mode == self.READ:
            return gzip.open(self.file_name, "rb")
        elif mode == self.WRITE:
            return gzip.open(self.file_name, "wb")
    
    def _get_bzip2(self):
        if mode == self.READ:
            return bz2.BZ2File(self.file_name, "r")
        elif mode == self.WRITE:
            return bz2.BZ2File(self.file_name, "w")

class FastqWriter(object):
    def __init__(self, file_name, compression=GUESS_COMPRESSION):
        self.fastq_io = FastqIO(
                                file_name,
                                mode=FastqIO.WRITE,
                                compression=compression
                               )
    
    def write_read(self, read):
        for line in read:
            self.fastq_io.write(line + "\n")

class FastqIterableReader(object):
    def __init__(self, file_name, mode=READ, compression=GUESS_COMPRESSION):
        self.fastq_io = FastqIO(file_name, mode, compression).get()
        self.num_records = 0
    
    def __iter__(self):
        return self
    
    def next(self):
        this_read = []
        for line in self.fastq_io:
            line = line.strip()
            if self.num_records < 0:
                if line[0] != "@":
                    continue
                else:
                    self.num_records += 1
            this_read.append(line)
            if len(this_read) == 4:
                return this_read
                self.num_records += 1
        raise StopIteration

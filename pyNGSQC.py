import gzip
import bz2
import os.path


#TODO
## Make proper docstrings

# Global option Variables
## Compression
GUESS_COMPRESSION = 0
NO_COMPRESSION = 1
GZIPPED = 2
BZIPP2ED = 3

## Barcode postion
BARCODE_FORWARD_ONLY = 0
BARCODE_REVERSE_ONLY = 0  # Not Implemented
BARCODE_FORWARD_REVERSE = 0  # Not Implemented

## CSV DIALECTS
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
                    "Invalid quality score %i from phred %s (ord %i)" %
                    (qual, phred, ord(phred))
                )
        return qual

    # General Methods:

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


class _GenericIO(object):
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
            raise ValueError(
                             "%i is not a valid compression mode" %
                              compression
                            )

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
            return self._get_plaintext()
        elif self.compression == GZIPPED:
            return self._get_gzip()
        elif self.compression == BZIPP2ED:
            return self._get_bzip2()
        else:
            raise ValueError(self.compression)

    def _get_plaintext(self):
        if self.mode == self.READ:
            return open(self.file_name, "rb")
        elif self.mode == self.WRITE:
            return open(self.file_name, "wb")

    def _get_gzip(self):
        if self.mode == self.READ:
            return gzip.open(self.file_name, "rb")
        elif self.mode == self.WRITE:
            return gzip.open(self.file_name, "wb")

    def _get_bzip2(self):
        if self.mode == self.READ:
            return bz2.BZ2File(self.file_name, "r")
        elif self.mode == self.WRITE:
            return bz2.BZ2File(self.file_name, "w")


class FastqIO():
    def __del__(self):
        self.io.close()


class FastqWriter(FastqIO):
    def __init__(self, file_name, compression=GUESS_COMPRESSION):
        self.io = _GenericIO(
                              file_name,
                              mode=_GenericIO.WRITE,
                              compression=compression
                            ).get()

    def write(self, read):
        for line in read:
            self.io.write(line + "\n")


class FastqReader(FastqIO):
    def __init__(
                  self,
                  file_name,
                  compression=GUESS_COMPRESSION
                ):
        self.io = _GenericIO(
                              file_name,
                              mode=_GenericIO.READ,
                              compression=compression
                            ).get()
        self.num_records = 0

    def __iter__(self):
        return self

    def next(self):
        this_read = []
        at_start = True
        for line in self.io:
            line = line.strip()
            if at_start:
                if line[0] != "@":
                    continue
                else:
                    at_start = False
            this_read.append(line)
            if len(this_read) == 4:
                return this_read
                self.num_records += 1
        raise StopIteration
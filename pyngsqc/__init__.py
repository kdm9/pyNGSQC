# Copyright 2012 Kevin Murray
# program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import gzip
import bz2
import os.path
from array import array

# TODO
## Make proper docstrings


# Global option Variables
## Defaults
DEFAULT_QUAL_OFFSET = 33
DEFAULT_QUAL_THRESHOLD = 15

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

## OUTPUT OBJECTS
STDOUT = 0

AMBIGUITY_DICT = {
    "A": ["A", "W", "M", "R", "N", "D", "H", "V"],
    "C": ["C", "S", "M", "Y", "N", "B", "H", "V"],
    "G": ["G", "S", "K", "R", "N", "B", "D", "V"],
    "T": ["T", "U", "W", "K", "Y", "N", "B", "D", "H"],
    "U": ["T", "U", "W", "K", "Y", "N", "B", "D", "H"],
    "N": ["A", "G", "C", "T", "U", "W", "K", "M", "R", "Y", "N", "B",
          "D", "H", "V"],
    "B": ["C", "G", "T", "U"],
    "D": ["A", "G", "T", "U"],
    "H": ["A", "C", "T", "U"],
    "V": ["A", "C", "G"],
    "K": ["G", "T", "U"],
    "M": ["A", "C"],
    "S": ["G", "C"],
    "W": ["A", "T", "U"],
    "Y": ["C", "T", "U"],
    "R": ["A", "G"]
}


class Base(object):

    num_good_reads = 0L
    num_bad_reads = 0L
    num_reads = 0L
    stats = {}

    def __init__(
            self,
            in_file_name,
            out_file_name,
            verbose=False,
            compression=GUESS_COMPRESSION,
            deduplicate_header=True,
            print_summary=False
    ):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        # self.compression = compression  # this is set in the reader and
        # writer
        self.verbose = verbose
        if self.in_file_name is not None:
            self.reader = FastqReader(
                self.in_file_name,
                compression=compression,
                deduplicate_header=deduplicate_header
            )
        if self.out_file_name is not None:
            self.writer = FastqWriter(
                self.out_file_name,
                compression=compression
            )
        self.print_summary = print_summary


class QualBase(Base):

    def __init__(
            self,
            # Inherited args
            in_file_name,
            out_file_name,
            # Local args
            # Inherited kwargs
            verbose=False,
            print_summary=False,
            compression=GUESS_COMPRESSION,
            deduplicate_header=True,
            # Local kwargs
            qual_offset=DEFAULT_QUAL_OFFSET,
            qual_threshold=DEFAULT_QUAL_THRESHOLD
    ):
        super(QualBase, self).__init__(
            in_file_name,
            out_file_name,
            verbose=verbose,
            print_summary=print_summary,
            compression=compression,
            deduplicate_header=True
        )
        self.qual_offset = qual_offset
        self.qual_threshold = qual_threshold


class _GenericFileHandle(object):
    READ = 0
    WRITE = 1

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


class _IOObject(object):

    def __init__(self):
        self.stats = {}
        self.stats["num_reads"] = 0

    def close(self):
        self.io.close()


class _Writer(_IOObject):

    def __init__(self, file_name, compression=GUESS_COMPRESSION):
        super(_Writer, self).__init__()
        self.io = _GenericFileHandle(
            file_name,
            mode=_GenericFileHandle.WRITE,
            compression=compression
        ).get()


class FastqWriter(_Writer):

    def __init__(self, file_name, compression=GUESS_COMPRESSION):
        super(FastqWriter, self).__init__(file_name, compression)

    def write(self, reads):
        if len(reads) == 4:
            self.stats["num_reads"] += 1
            for line in reads:
                self.io.write(line + "\n")
        elif len(reads) % 4 == 0:
            self.stats["num_reads"] += len(reads) / 4
            for line in reads:
                self.io.write(line + "\n")
        else:
            raise ValueError("Bad Fastq Read: %s" % repr(reads))


class FastaWriter(_Writer):
    def __init__(self, file_name, compression=GUESS_COMPRESSION):
        super(FastaWriter, self).__init__(file_name, compression)

    def write(self, reads):
        if len(reads) == 2:
            self.stats["num_reads"] += 1
            for line in reads:
                self.io.write(line + "\n")
        elif len(reads) % 2 == 0:
            self.stats["num_reads"] += len(reads) / 2
            for line in reads:
                self.io.write(line + "\n")
        else:
            raise ValueError("Bad Fasta Read: %s" % repr(reads))


class _Reader(_IOObject):

    def __init__(
            self,
            file_name,
            deduplicate_header=True,
            compression=GUESS_COMPRESSION
    ):
        super(_Reader, self).__init__()
        self.file_name = file_name
        self.deduplicate_header = deduplicate_header
        self.io = _GenericFileHandle(
            self.file_name,
            mode=_GenericFileHandle.READ,
            compression=compression
        ).get()
        self.stats["num_reads"] = 0L

    def __iter__(self):
        return self


class FastqReader(_Reader):

    def __init__(
            self,
            file_name,
            deduplicate_header=True,
            compression=GUESS_COMPRESSION
    ):
        super(FastqReader, self).__init__(
            file_name,
            deduplicate_header,
            compression
        )

    def next(self):
        this_read = []
        at_start = True
        num_lines = 0
        for line in self.io:
            num_lines += 1
            line = line.strip()
            if at_start:
                if line[0] != "@":
                    continue
                else:
                    at_start = False
            this_read.append(line)
            if len(this_read) == 4:
                if not len(this_read[1]) == len(this_read[3]):
                    err = "Read %s has seq and qual of different lengths"
                    raise ValueError(err % repr(this_read))
                if not this_read[2][0] == "+":
                    err = "Read %s has no quality header, or is misformed"
                    raise ValueError(err % repr(this_read))
                if self.deduplicate_header:
                    # Save disk space, remove duplicate headers
                    if this_read[0][1:] == this_read[2][1:]:
                        this_read[2] = "+"
                self.stats["num_reads"] += 1
                return this_read
        raise StopIteration


class FastqRandomAccess(_Reader):

    def __init__(
            self,
            file_name,
            deduplicate_header=True,
            compression=GUESS_COMPRESSION
    ):
        super(FastqRandomAccess, self).__init__(
            file_name,
            deduplicate_header,
            compression
        )
        self.record_positions = array("L")
        self._build_cache()

    def get(self, index):
        start_pos = self.record_positions[index]
        self.io.seek(start_pos)
        this_read = []
        for iii in xrange(4):
            this_read.append(self.io.readline())
        # print index, this_read
        if not len(this_read[1]) == len(this_read[3]):
            err = "Read %s has seq and qual of different lengths"
            raise ValueError(err % repr(this_read))
        if not this_read[2][0] == "+":
            err = "Read %s has no quality header, or is misformed"
            raise ValueError(err % repr(this_read))
        if self.deduplicate_header:
            # Save space, remove duplicate headers
            if this_read[0][1:] == this_read[0][2:]:
                this_read[2] = "+"
        self.stats["num_reads"] += 1
        return this_read

    def _build_cache(self):
        at_start = True
        current_pos = 0
        record_str = ""
        for line in self.io:
            if at_start:
                if line[0] != "@":
                    current_pos += len(line)
                    continue
                else:
                    at_start = False
            if line[0] == "@":
                current_pos += len(record_str)
                self.record_positions.append(current_pos)
                record_str = line
            else:
                record_str += line


def base_match(base_1, base_2, allow_ambiguity=True):
    """
    base_match(base_1, base_2):
        Returns true if base_1 and base_2 are equivalent according to IUPAC
        rules
    """
    if base_1 == base_2:
        return True
    if allow_ambiguity:
        try:
            result = base_2 in AMBIGUITY_DICT[base_1]
        except KeyError:
            raise ValueError("(%s, %s) is an invalid base pair" % (
                base_1,
                base_2
            )
            )
    else:
        result = base_2 == base_1
    return result


def seq_match(seq_1, seq_2, mismatches=0, allow_ambiguity=False):
    if seq_1 == seq_2:
        return True
    if mismatches < 1:
        return False
    if allow_ambiguity:
        this_mismatches = 0
        for iii in xrange(len(seq_1)):
            if not base_match(seq_1[iii], seq_2[iii], allow_ambiguity=True):
                this_mismatches += 1
            if this_mismatches > mismatches:
                return False
        return this_mismatches <= mismatches
    else:
        for iii in xrange(len(seq_1)):
            if seq_1[iii] != seq_2[iii]:
                mismatches -= 1
                if mismatches < 0:
                    return False
        return True


def num_Ns_in_read(read):
    seq = read[1]
    seq = seq.upper()
    if "N" not in seq:
        return False
    else:
        return seq.count('N')


def get_qual_from_phred(phred, offset):
    qual = ord(phred) - offset
    if qual < 0:
        raise ValueError(
            "Invalid quality score %i from phred %s (ord %i)" %
            (qual, phred, ord(phred))
        )
    return qual


def get_phred_from_qual(qual, offset):
    return chr(qual + offset)


def percentile_from_counts(count_list, percentile):
    """Returns the median value from a list whose values are counts of the
    index i
    """
    # The index of the median, or of the lower of the two values if sum is even
    halfway_index = int(round(float(sum(count_list) * percentile))) - 1
    current_index = 0
    pos_within_value = 0  # Governs when we skip to the next median value
    median = 0
    while current_index <= halfway_index:
        # If we have not iterated through all counts of this value of median
        if pos_within_value < count_list[median]:
            pos_within_value += 1
            current_index += 1
        else:
            # Move to next median
            median += 1
            pos_within_value = 0
    # At this point, median == lower of two values if the number of counts is
    # even, so we need to average the current value of median with the
    # following value of median (which may be the same number) to get the
    # median. If the number of counts is odd, then median == true median
    if sum(count_list) % 2 == 0:
        lower_score = median
        if count_list[median] < pos_within_value + 1:
            upper_score = median + 1
        else:
            upper_score = median
        median = float(lower_score + upper_score) / 2.0
    return float(median)


def whiskers_from_counts(count_list):
    q1 = percentile_from_counts(count_list, 0.25)
    q3 = percentile_from_counts(count_list, 0.75)
    median = percentile_from_counts(count_list, 0.5)
    iqr = q3 - q1
    left_whisker = median - 1.5 * iqr
    right_whisker = median + 1.5 * iqr
    return (left_whisker, right_whisker)


def dict_to_tuples(this_dict):
    tuples = []
    for key, value in this_dict:
        tuples.append((key, value))
    return tuples


def tuples_to_dict(tuples):
    this_dict = {}
    for key, value in tuples:
        this_dict[key] = value
    return this_dict


def convert_phred_offset(in_phred, in_qual_offset, out_qual_offset):
    out_phred = ""
    for char in in_phred:
        out_phred += get_phred_from_qual(
            get_qual_from_phred(char, in_qual_offset),
            out_qual_offset
        )
    return out_phred

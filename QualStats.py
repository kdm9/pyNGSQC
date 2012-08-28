#from sys import stderr
import pyNGSQC


class QualStats(pyNGSQC.NGSQC):
    def __init__(
                 self,
                 in_file_name,
                 qual_offset=64,
                 min_length=15,
                 output=pyNGSQC.STDOUT,
                 compression=pyNGSQC.GUESS_COMPRESSION
                ):
        self.in_file_name = in_file_name
        self.output = output
        self.reader = pyNGSQC.FastqReader(
                                           self.in_file_name,
                                           compression=compression
                                         )
        self.qual_offset = qual_offset
        self.min_length = min_length
        self.num_reads = 0
        self.position_bases = []
        self.position_scores = []

        # Set default dicts for the position_bases and position_qualities lists
        bases = list("AGCTSWRYKAN-")
        self.initial_base_counts = dict.fromkeys(bases, 0)
        # 126 is the last printable character in ASCII, store as list as keys
        # are integers, saves mem
        self.initial_score_counts = [0 for iii in xrange(self.qual_offset, 126)]

    def print_summary(self):
        import json
        print json.dumps(self.position_scores, indent=2)
        print json.dumps(self.position_bases, indent=2)

    def process_read(self, read):
        while len(self.position_bases) < 101:  # len(read[1]):
            self.position_bases.append(dict(self.initial_base_counts))
        while len(self.position_scores) < 101:  # len(read[1]):
            self.position_scores.append(list(self.initial_score_counts))

        for bbb in xrange(len(read[1])):
            base = read[1][bbb]
            score = self._get_qual_from_phred(read[3][bbb])
            self.position_scores[bbb][score] += 1
            self.position_bases[bbb][base] += 1

    def run(self):
        for read in self.reader:
            self.num_reads += 1
            self.process_read(read)
        self.print_summary()
        return True

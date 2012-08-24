from sys import stderr
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
        # 126 is the last printable character in ASCII
        scores = [iii for iii in xrange(self.qual_offset, 126)]
        self.initial_score_counts = dict.fromkeys(scores, 0)

    def print_output(self):
        stderr.write("QualTrimmer finished:\n")
        stderr.write(
                      "\t%i sequences passed QC, wrote them to %s\n" %
                      (self.num_good_reads, self.out_file_name)
                    )
        stderr.write(
                      "\t%i sequences failed QC, and were ignored\n" %
                      self.num_bad_reads,
                    )

    def process_read(self, read):
        while len(self.position_scores) < len(read[1]):
            self.position_bases.append(dict(self.initial_base_counts))
        while len(self.position_scores) < len(read[1]):
            self.position_scores.append(dict(self.initial_score_counts))

        for bbb in xrange(len(read[1])):
            base = read[1][bbb]
            score = self._get_qual_from_phred(read[3][bbb])
            try:
                self.position_scores[bbb][score] += 1
                self.position_bases[bbb][base] += 1
            except KeyError:
                pass

    def run(self):
        for read in self.reader:
            self.num_reads += 1
            self.process_read(read)
        import json
        print json.dumps(self.position_scores, indent=2)
        print json.dumps(self.position_bases, indent=2)
        self.print_summary()
        return True


qs = QualStats("/home/kevin/UniWork/BIOL3157/Assignments/2/U4852380.txt")
qs.run()

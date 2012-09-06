#from sys import stderr
import pyNGSQC
from copy import deepcopy


class _ReadTask():

    def __init__(self):
        super(_ReadTask, self).__init__()


class QualStats(pyNGSQC.NGSQC):
    def __init__(
                 self,
                 in_file_name,
                 qual_offset=64,
                 min_length=15,
                 verbose=False,
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
        self.positions = []

        # Set default dicts for the position_bases and position_qualities lists
        bases = list("AGCTN")
        self.columns = [
            "count",
            "min",
            "max",
            "sum",
            "mean",
            "Q1",
            "median",
            "Q3",
            "IQR",
            "lW",
            "rW",
            "A",
            "C",
            "G",
            "T",
            "N",
            "total",
            "GC"
            ]
        self.initial_dict = {
            "bases": dict.fromkeys(bases, 0),
            # 126 is the last printable character in ASCII, store as list as
            # keys are integers, saves mem
            "scores": [0 for iii in xrange(self.qual_offset, 126)],
            # Dict of summary data to be printed
            "summary": dict.fromkeys(self.columns, 0.0)
            }

    def _print_summary(self):
        self._summarize_data()
        print "\t".join(self.columns)
        for position in self.positions:
            values = [str(position["summary"][col]) for col in self.columns]
            print "\t".join(values)

    def _summarize_data(self):
        for position in self.positions:
            # First, calculate summary stats on scores
            position["summary"]["min"] = float(min(position["scores"]))
            position["summary"]["max"] = float(max(position["scores"]))
            # Sum of the list is the total number of scores
            position["summary"]["count"] = float(sum(position["scores"]))
            # The actual sum is the sum of (index * count)
            for sss in xrange(len(position["scores"])):
                position["summary"]["sum"] += \
                 float(sss * position["scores"][sss])
            # Mean
            position["summary"]["mean"] = position["summary"]["sum"] / \
             position["summary"]["count"]
            # Quartiles, median, iqr and whiskers
            position["summary"]["Q1"] = pyNGSQC._percentile_from_counts(
                position["scores"],
                0.25
                )
            position["summary"]["median"] = pyNGSQC._percentile_from_counts(
                position["scores"],
                0.5
                )
            position["summary"]["Q3"] = pyNGSQC._percentile_from_counts(
                position["scores"],
                0.75
                )
            position["summary"]["IQR"] = position["summary"]["Q3"] - \
             position["summary"]["Q1"]

            (position["summary"]["lW"], position["summary"]["rW"]) = \
             pyNGSQC._whiskers_from_counts(position["scores"])

            # Calculate Base Counts
            total_count = 0
            for base, count in position["bases"].items():
                total_count += count
                position["summary"][base] = count
            g_plus_c = position["summary"]["G"] + position["summary"]["C"]
            position["summary"]["total"] = total_count
            position["summary"]["GC"] = float(g_plus_c) / float(total_count)

    def _process_read(self, read):
        while len(self.positions) < len(read[1]):
            self.positions.append(deepcopy(self.initial_dict))

        for pos in xrange(len(read[1])):
            base = read[1][pos]
            score = self._get_qual_from_phred(read[3][pos])
            self.positions[pos]["scores"][score] += 1
            try:
                self.positions[pos]["bases"][base] += 1
            except KeyError:
                self.positions[pos]["bases"]["N"] += 1

    def run(self):
        for read in self.reader:
            self.num_reads += 1
            self._process_read(read)
        self._print_summary()
        return True

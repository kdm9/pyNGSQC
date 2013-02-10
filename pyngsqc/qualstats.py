#Copyright 2012 Kevin Murray
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pyngsqc
from copy import deepcopy


class _QualStatsTask():

    def __init__(self):
        super(_QualStatsTask, self).__init__()


class QualStats(pyngsqc.QualBase):

    def __init__(
            self,
            # Inherited args
            in_file_name,
            # Local args
            # Inherited kwargs
            qual_offset=pyngsqc.DEFAULT_QUAL_OFFSET,
            compression=pyngsqc.GUESS_COMPRESSION,
            deduplicate_header=True,
            verbose=False,
            print_summary=False,
            # Local kwargs
            output=pyngsqc.STDOUT
            ):
        # Initialise base class
        super(QualStats, self).__init__(
                in_file_name,
                out_file_name = None,  # This class doesn't have one
                qual_offset=qual_offset,
                compression=compression,
                deduplicate_header=deduplicate_header,
                verbose=verbose,
                print_summary=print_summary
                )
        # Initialise local variables
        self.output = output
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
                # 126 is the last printable character in ASCII, store as list
                # as keys are integers, saves mem
                "scores": [0 for iii in xrange(self.qual_offset, 126)],
                # Dict of summary data to be printed
                "summary": dict.fromkeys(self.columns, 0.0)
                }

    def _print_summary(self):
        print
        print "cycle\t", "\t".join(self.columns)
        cycle = 0
        for position in self.positions:
            cycle += 1
            values = [str(position["summary"][col]) for col in self.columns]
            print "%i\t" % cycle, "\t".join(values)

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
            position["summary"]["Q1"] = pyngsqc.percentile_from_counts(
                    position["scores"],
                    0.25
                    )
            position["summary"]["median"] = pyngsqc.percentile_from_counts(
                    position["scores"],
                    0.5
                    )
            position["summary"]["Q3"] = pyngsqc.percentile_from_counts(
                    position["scores"],
                    0.75
                    )
            position["summary"]["IQR"] = position["summary"]["Q3"] - \
             position["summary"]["Q1"]

            (position["summary"]["lW"], position["summary"]["rW"]) = \
             pyngsqc.whiskers_from_counts(position["scores"])

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
            score = pyngsqc.get_qual_from_phred(read[3][pos], self.qual_offset)
            self.positions[pos]["scores"][score] += 1
            try:
                self.positions[pos]["bases"][base] += 1
            except KeyError:
                self.positions[pos]["bases"]["N"] += 1

    def run(self):
        for read in self.reader:
            self.num_reads += 1
            self._process_read(read)
        self._summarize_data()
        if self.print_summary:
            self._print_summary()
        return True

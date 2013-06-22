# Copyright 2012 Kevin Murray
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from sys import stderr
import pyngsqc
import csv
from os import path
from pyngsqc import seq_match, _parallel

FORWARD_ONLY = 0
REVERSE_ONLY = 1
FORWARD_OR_REVERSE = 2


class BarcodeSplitter(pyngsqc.Base):
    def __init__(
            self,
            # Inherited args
            in_file_name,
            output_dir,
            # Local args
            barcode_file,
            # Inherited kwargs
            compression="guess",
            deduplicate_header=True,
            verbose=False,
            print_summary=False,
            # Local kwargs
            barcode_end=FORWARD_ONLY,
            mismatches=0,
            write_to_header=False
    ):
        # Initialise base class
        super(BarcodeSplitter, self).__init__(
            in_file_name,
            out_file_name=None,  # we dont have one for this class
            compression=compression,
            deduplicate_header=deduplicate_header,
            verbose=verbose,
            print_summary=print_summary
        )
        self.output_dir = output_dir
        self.mismatches = mismatches
        self.barcodes = self._set_barcodes_from_file(
            barcode_file,
            pyngsqc.SNIFF_CSV_DIALECT
        )
        self.write_to_header = write_to_header
        self.writer = _BarcodeWriter(
            self.barcodes,
            self.in_file_name,
            self.output_dir
        )

    def _sniff_csv_dialect(self, file_name):
        csv_fh = open(file_name, "rb")
        sample = ""
        for i in range(10):
            sample += bytes(csv_fh.readline()).decode("UTF-8")
        csv_dialect = csv.Sniffer().sniff(sample)
        csv_fh.close()
        return csv_dialect

    def _set_barcodes_from_file(self, file_name, csv_dialect=csv.excel):
        barcodes = {}
        if csv_dialect == pyngsqc.SNIFF_CSV_DIALECT:
            csv_dialect = self._sniff_csv_dialect(file_name)
        csv_fh = open(file_name)
        csv_reader = csv.reader(csv_fh, dialect=csv_dialect)
        for line in csv_reader:
            barcodes[line[0]] = "".join(line[1:])
        csv_fh.close()
        return barcodes

    def _parse_read_barcode(self, read):
        for barcode in self.barcodes:
            barcode_len = len(barcode)
            read_barcode = read[1][0:barcode_len]
            if seq_match(barcode, read_barcode, mismatches=self.mismatches):
                if self.write_to_header:
                    read[0] += " bcd:%s desc:%s" % (
                        barcode,
                        self.barcodes[barcode]
                    )
                read[1] = read[1][barcode_len:]
                read[3] = read[3][barcode_len:]
                self.writer.write((barcode, read))

    def _print_summary(self):
        stderr.write("Barcode Splitter finished:\n")
        stderr.write(
            "\t%i sequences analysed, containing %s\n" % (
            self.stats["writer"]["num_reads"],
            len(self.stats["writer"]["barcode_counts"])
            )
        )

        if self.verbose:
            stderr.write("\tThe following barcodes were parsed:\n")
            for bcd, count in self.stats["writer"]["barcode_counts"].items():
                stderr.write("\t%s:\t%i\n" % (bcd, count))

    def run(self):
        if len(self.barcodes) < 1:
            raise ValueError(
                "You must supply a barcode dict or file before"
                " run()-ing BarcodeSplitter"
            )
        for read in self.reader:
            self._parse_read_barcode(read)

        self.stats["reader"] = self.reader.stats
        self.stats["writer"] = self.writer.stats
        if self.print_summary:
            self._print_summary()

    def run_parallel(self):
        if len(self.barcodes) < 1:
            raise ValueError(
                "You must supply a barcode dict or file before"
                " run()-ing BarcodeSplitter"
            )
        runner = _parallel.ParallelRunner(
            BarcodeSplitTask,
            self.reader,
            self.writer,
            (  # Task options
                self.barcodes,
                self.mismatches,
                self.write_to_header,
            )
        )
        runner.run()

        self.stats["runner"] = runner.stats
        self.stats["reader"] = self.reader.stats
        self.stats["writer"] = self.writer.stats
        if self.print_summary:
            self._print_summary()


class BarcodeSplitTask(object):

    def __init__(self, read, barcodes, mismatches, write_to_header):
        self.read = read
        self.barcodes = barcodes
        self.mismatches = mismatches
        self.write_to_header = write_to_header

    def __call__(self):
        for barcode in self.barcodes:
            barcode_len = len(barcode)
            read_barcode = self.read[1][0:barcode_len]
            if seq_match(barcode, read_barcode, mismatches=self.mismatches):
                if self.write_to_header:
                    self.read[0] += " bcd:%s desc:%s" % (
                        barcode,
                        str(self.barcodes[barcode])
                    )
                self.read[1] = self.read[1][barcode_len:]
                self.read[3] = self.read[3][barcode_len:]
                return (barcode, self.read)
        return (None, self.read)


class _BarcodeWriter(pyngsqc.Base):
    """Provides a pyngsqc.FastqWriter compatible interface to write to many
    files at once.
    """
    def __init__(self, barcodes, in_file_name, output_dir=None):
        self.in_file_name = in_file_name
        self.output_dir = output_dir
        self.barcodes = barcodes
        self.stats = {}
        self.stats["barcode_counts"] = {}
        self.barcode_files = {}
        self.stats["num_reads"] = 0

    def _get_barcode_writer(self, barcode):
        if self.barcodes[barcode] == "" or self.barcodes[barcode] is None:
            identifier = barcode
        else:
            # Not sure if this is a good idea, is possibility of the
            # barcode's decsription being non-unique
            identifier = self.barcodes[barcode]
        if self.output_dir is None:
            dir_path = path.dirname(self.in_file_name)
        else:
            dir_path = self.output_dir
        split_path = path.basename(self.in_file_name).split(".")
        out_path = path.join(dir_path, split_path[0] + "_%s" % identifier)
        if len(split_path) > 1:
            # If the path had extensions, add them
            out_path += "." + ".".join(split_path[1:])
        writer = pyngsqc.FastqWriter(out_path)
        return writer

    def write(self, pair):
        barcode, read = pair
        if barcode is not None:
            if barcode not in self.barcode_files:
                self.barcode_files[barcode] = \
                    self._get_barcode_writer(barcode)
            self.stats["num_reads"] += 1
            self.barcode_files[barcode].write(read)
            try:
                self.stats["barcode_counts"][barcode] += 1
            except KeyError:
                self.stats["barcode_counts"][barcode] = 1

    def close(self):
        for barcode in self.barcode_files:
            self.barcode_files[barcode].close()


class PairedBarcodeSplitter(BarcodeSplitter):
    def __init__(
            self,
            pair_1_file_name,
            pair_2_file_name,
            output_dir=None,
            compression="guess",
            verbose=False
    ):
        self.pair_1_file_name = pair_1_file_name
        self.pair_2_file_name = pair_2_file_name
        self.output_dir = output_dir
        self.pair_1_reader = pyngsqc.FastqReader(
            self.pair_1_file_name,
            compression=compression
        )
        self.pair_2_reader = pyngsqc.FastqReader(
            self.pair_2_file_name,
            compression=compression
        )
        self.verbose = verbose
        self.barcodes = {}
        self.barcode_counts = {}
        # Note that files will be stored as a tuple of writers (pair_1, pair_2)
        self.barcode_files = {}
        self.num_reads = 0

    def _get_barcode_writers(self, barcode):
        if self.barcodes[barcode] == "" or self.barcodes[barcode] is None:
            identifier = barcode
        else:
            # Not sure if this is a good idea, is possibility of the
            # barcode's decsription being non-unique
            identifier = self.barcodes[barcode]
        if self.output_dir is None:
            dir_path = path.dirname(self.pair_1_file_name)
        else:
            dir_path = self.output_dir

        # Read 1
        read_1_split_path = path.basename(self.pair_1_file_name).split(".")
        out_path = path.join(
            dir_path,
            read_1_split_path[0] + "_%s" % identifier
        )
        if len(read_1_split_path) > 1:
            # If the path had extensions, add them
            out_path += "." + ".".join(read_1_split_path[1:])
        read_1_writer = pyngsqc.FastqWriter(out_path)

        # Read 2
        read_2_split_path = path.basename(self.pair_1_file_name).split(".")
        out_path = path.join(
            dir_path,
            read_2_split_path[0] + "_%s" % identifier
        )

        # If the path had extensions, add them
        if len(read_2_split_path) > 1:
            out_path += "." + ".".join(read_2_split_path[1:])
        read_2_writer = pyngsqc.FastqWriter(out_path)

        return (read_1_writer, read_2_writer)

    def _parse_paired_reads_barcode(self, paired_reads):
        (read_1, read_2) = paired_reads
        for barcode in self.barcodes:
            barcode_len = len(barcode)
            if read_1[1][0:barcode_len] == barcode:
                if barcode not in self.barcode_files:
                    writers = self._get_barcode_writers(barcode)
                    self.barcode_files[barcode] = writers
                if barcode not in self.barcode_counts:
                    self.barcode_counts[barcode] = 0
                read_1[0] += " bcd:%s desc:%s" % \
                    (barcode, self.barcodes[barcode])
                read_1[1] = read_1[1][barcode_len:]
                read_1[3] = read_1[3][barcode_len:]
                (read_1_writer, read_2_writer) = self.barcode_files[barcode]
                read_1_writer.write(read_1)
                read_2_writer.write(read_2)
                self.barcode_counts[barcode] += 1

    def run(self):
        if len(self.barcodes) < 1:
            raise ValueError(
                "You must supply a barcode dict or file before"
                " run()-ing BarcodeSplitter"
            )
        for paired_reads in zip(self.pair_1_reader, self.pair_2_reader):
            self.num_reads += 1
            self._parse_paired_reads_barcode(paired_reads)

        self.stats["reader"] = self.reader.stats
        self.stats["writer"] = self.writer.stats
        if self.print_summary:
            self._print_summary()

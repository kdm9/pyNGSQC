from sys import stderr
import pyNGSQC
import csv
from os import path
from pyNGSQC import seq_match
FORWARD_ONLY = 0
REVERSE_ONLY = 1
FORWARD_OR_REVERSE = 2


class BarcodeSplitter(pyNGSQC.NGSQC):
    def __init__(
            self,
            in_file_name,
            output_dir=None,
            barcode_end=FORWARD_ONLY,
            compression=pyNGSQC.GUESS_COMPRESSION,
            mismatch=1,
            verbose=False
            ):
        self.in_file_name = in_file_name
        self.output_dir = output_dir
        self.reader = pyNGSQC.FastqReader(
            self.in_file_name,
            compression=compression
            )
        self.mismatch = mismatch
        self.verbose = verbose
        self.barcodes = {}
        self.barcode_counts = {}
        self.barcode_files = {}
        self.num_reads = 0

    def _sniff_csv_dialect(self, file_name):
        csv_fh = open(file_name, "rb")
        sample = ""
        for i in xrange(10):
            sample += csv_fh.readline()
        csv_dialect = csv.Sniffer().sniff(sample)
        csv_fh.close()
        return csv_dialect

    def set_barcodes_from_file(self, file_name, csv_dialect=csv.excel):
        """
        """
        if csv_dialect == pyNGSQC.SNIFF_CSV_DIALECT:
            csv_dialect = self._sniff_csv_dialect(file_name)
        csv_fh = open(file_name, "rb")
        csv_reader = csv.reader(csv_fh, dialect=csv_dialect)
        for line in csv_reader:
            self.barcodes[line[0]] = "".join(line[1:])
        #csv_reader.close()
        csv_fh.close()

    def set_barcodes_from_dict(self, barcodes):
        """
        """
        self.barcodes = barcodes

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
        writer = pyNGSQC.FastqWriter(out_path)
        return writer

    def _parse_read_barcode(self, read):
        for barcode in self.barcodes:
            barcode_len = len(barcode)
            read_barcode = read[1][0:barcode_len]
            if seq_match(barcode, read_barcode, mismatches=1):
            #if barcode == read_barcode:
                if barcode not in self.barcode_files:
                    writer = self._get_barcode_writer(barcode)
                    self.barcode_files[barcode] = writer
                if barcode not in self.barcode_counts:
                    self.barcode_counts[barcode] = 0
                read[0] += " bcd:%s desc:%s" % (
                    barcode,
                    self.barcodes[barcode]
                    )
                read[1] = read[1][barcode_len:]
                read[3] = read[3][barcode_len:]
                self.barcode_files[barcode].write(read)
                self.barcode_counts[barcode] += 1

    def print_summary(self):
        stderr.write("Barcode Splitter finished:\n")
        stderr.write(
            "\t%i sequences analysed, containing %s\n" %
            (self.num_reads, len(self.barcode_counts))
            )

        if self.verbose:
            stderr.write("\tThe following barcodes were parsed:\n")
            for barcode, count in self.barcode_counts.iteritems():
                stderr.write("\t%s:\t%i\n" % (barcode, count))

    def run(self):
        if len(self.barcodes) < 1:
            raise ValueError(
                "You must supply a barcode dict or file before" +
                " run()-ing BarcodeSplitter"
                )
            return False
        for read in self.reader:
            self.num_reads += 1
            self._parse_read_barcode(read)
        self.print_summary()
        return True


class PairedBarcodeSplitter(BarcodeSplitter):
    def __init__(
                 self,
                 pair_1_file_name,
                 pair_2_file_name,
                 output_dir=None,
                 compression=pyNGSQC.GUESS_COMPRESSION,
                 verbose=False
                ):
        self.pair_1_file_name = pair_1_file_name
        self.pair_2_file_name = pair_2_file_name
        self.output_dir = output_dir
        self.pair_1_reader = pyNGSQC.FastqReader(
                                                 self.pair_1_file_name,
                                                 compression=compression
                                                )
        self.pair_2_reader = pyNGSQC.FastqReader(
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
        read_1_writer = pyNGSQC.FastqWriter(out_path)

        # Read 2
        read_2_split_path = path.basename(self.pair_1_file_name).split(".")
        out_path = path.join(
            dir_path,
            read_2_split_path[0] + "_%s" % identifier
            )

        # If the path had extensions, add them
        if len(read_2_split_path) > 1:
            out_path += "." + ".".join(read_2_split_path[1:])
        read_2_writer = pyNGSQC.FastqWriter(out_path)

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
                "You must supply a barcode dict or file before" +
                " run()-ing BarcodeSplitter"
                )
        for paired_reads in zip(self.pair_1_reader, self.pair_2_reader):
            self.num_reads += 1
            self._parse_paired_reads_barcode(paired_reads)
        self.print_summary()
        return True

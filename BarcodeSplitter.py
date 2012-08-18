from sys import stderr 
import gzip
import bz2
import pyNGSQC
import csv

class BarcodeSplitter(pyNGSQC.NGSQC):
    def __init__(
                 self,
                 in_file_name,
                 out_file_name,
                 barcode_end=pyNGSQC.BARCODE_FORWARD_ONLY, 
                 append=False,
                 remove_trailing_Ns=True,
                 compression=pyNGSQC.GUESS_COMPRESSION
                ):
        
        super(BarcodeSplitter, self).__init__(
                                            in_file_name,
                                            out_file_name,
                                            qual_offset=qual_offset,
                                            append=append,
                                            compression=compression
                                           )
        self.barcode_map = {}
        self.barcode_counts = {}
        self.barcode_output_files = {}
    
    def set_barcodes_from_file(self, file_name, csv_dialect=csv.excel):
        csv_fh = open(file_name, "rb")
        if csv_dialect == pyNGSQC.SNIFF_CSV_DIALECT:
            sample = ""
            for i in xrange(10):
                sample += csv_fh.readline()
            csv_dialect = csv.Sniffer.sniff(sample)
            csv_fh.close()
            csv_fh = open(file_name, "rb")
        csv_reader = csv.reader(csv_fh, dialect=csv_dialect)
        for line in csv_fh:
            self.barcode_map = {
                                "ATCG":"Description"
                               }
        
    
    def set_barcodes_from_dict(self, dict):
        pass
    
    def pass_read_barcode(self, read):
        for barcode in self.barcode_map.keys():
            # 
            barcode_len = len(barcode)
            if read[1][0:barcode_len] == barcode:
                
                return self.barcode_map[barcode]
            pass
        
    def process_read(self, read):
        pass
    

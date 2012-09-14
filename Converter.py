# -*- coding: utf-8 *-*
import pyNGSQC


class FastqToFasta(pyNGSQC.Base):

    def __init__(
             self,
             in_file_name,
             out_file_name,
             remove_header=False,
             compression=pyNGSQC.GUESS_COMPRESSION
            ):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.reader = pyNGSQC.FastqReader(
            self.in_file_name,
            compression=compression
            )
        self.writer = pyNGSQC.FastaWriter(
            self.out_file_name,
            compression=compression
            )
        self.remove_header = remove_header

    def run(self):
        for read in self.reader:
            fasta_read = []
            if self.remove_header:
                header = ">%i" % self.reader.num_reads
            else:
                header = ">%s" % read[0][1:]  # Keep fastq header

            fasta_read.append(header)
            fasta_read.append(read[1])  # Seq

            self.writer.write(fasta_read)
        return True

#!/usr/bin/env python

"""libcharm.py: Provides codon harmonization functions to front-ends"""

from sys import exit

from urllib.request import Request, urlopen
from urllib.error import URLError

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from bs4 import BeautifulSoup


__version__ = "0.1"


def open_input_file(input):
    try:
        parser = SeqIO.parse(input, 'fasta', IUPAC.unambiguous_dna)
    except:
        try:
            parser = SeqIO.parse(input, 'fasta', IUPAC.unambiguous_rna)
        except:
            exit(1)
    content = []
    for record in parser:
        content.append(record)
    if len(content) > 1:
        print('More than one sequence found in file! Only the first sequence will be used!')

    seq = content[0].seq
    return seq


class CodonUsageTable():
    def __init__(self, url):
        self.url = url
        self.usage_table = {}
        self.fetch_codon_usage_table()


    def add_to_table(self, codon, aa, frequency):
        """Add codon and usage frequency to table"""

        if aa in self.usage_table:
        # If the aa is already present in the table, just add the new codon
            self.usage_table[aa][codon] = {'f': frequency}
        else:
            # Else, create a new entry for the aa and add the codon
            self.usage_table[aa] = {}
            self.usage_table[aa][codon] = {'f': frequency}

    def fetch_codon_usage_table(self):
        """Fetch the codon table from http://www.kazusa.or.jp/codon"""

        request = Request(self.url)
        try:
            opener = urlopen(request)
        except URLError as e:
            if hasattr(e, 'reason'):
                print('Failed to reach server: %s' % e.reason)
            elif hasattr(e, 'code'):
                print('Server responded with HTTP error code: %s' % e.code)
            else:
                print(opener)

        response = opener.read()
        soup = BeautifulSoup(response)

        table_string = str(soup.pre) # Parse the HTML response and look for <pre></pre>
        # section, containing the usage table

        table_string = table_string.replace('<pre>\n', '')  # remove <pre></pre> tags
        table_string = table_string.replace('\n</pre>', '') #

        table_lines = table_string.split('\n') # Split in lines
        for line in table_lines: # Iterate over the lines
            lines = line.split(')') # Splitting the lines at ")" will result in substrings representing the
            # a single codon each
            for codon_raw in lines:
                codon_raw = codon_raw.strip() # strip whitespace characters from the substring

                codon = codon_raw[:3].strip().replace('U', 'T') # The first three characters are the codon
                if codon:
                    aa = codon_raw[4:5].strip() # Position 5 is the aa in one letter code
                    frequency = float(codon_raw[6:10].strip()) # position 6 to 10 is the usage frequency;
                    # convert to float
                    self.add_to_table(codon, aa, frequency)


class Sequence():
    """Provides methods for storage and manipulation of sequences"""

    def __init__(self, sequence, origin_id, host_id):
        # Set translation table for original sequence
        self.translation_table = CodonTable.unambiguous_dna_by_name["Standard"]
        # Reformat and sanitize sequence string (remove whitespaces, change to uppercase)
        if type(sequence) is 'str':
            if 'U' in sequence:
                sequence = sequence.replace('U', 'T')
            self.original_sequence = Seq(''.join(sequence.upper().split()), IUPAC.unambiguous_dna)
        else:
            self.original_sequence = sequence
        self.original_translated_sequence = self.translate_sequence(self.original_sequence, cds=True)
        self.harmonized_sequence = ''
        self.codons = []
        self.split_to_codons()

        self.usage_origin = CodonUsageTable(
            'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={}&aa=1&style=N'.format(origin_id))
        self.usage_host = CodonUsageTable(
            'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={}&aa=1&style=N'.format(host_id))

        self.harmonize_codons(self.usage_origin, self.usage_host)
        self.construct_new_sequence()

    def translate_sequence(self, sequence, cds=True):
        try:
            translated_sequence = sequence.translate(table=self.translation_table, cds=cds)
        except:
            print("Sequence is not a valid CDS!")
            exit(1)
        return translated_sequence

    def get_harmonized_codons(self):
        harmonized_codons = []
        for codon in self.codons:
            if str(codon['original']) != str(codon['new']):
                harmonized_codons.append(codon)
        return harmonized_codons


    def split_to_codons(self):
        """Split the sequence into codons."""

        for codon in chunks(self.original_sequence, 3):
            self.codons.append({'original': codon,
                                'new': None,
                                'origin_f': None,
                                'target_f': None,
                                'initial_df': None,
                                'final_df': None,
                                'aa': str(codon.translate(table=self.translation_table))})

    def harmonize_codons(self, usage_origin, usage_target):
        for codon in self.codons:
            aa = codon['aa']
            orig_codon = str(codon['original'])

            origin_f = usage_origin.usage_table[aa][orig_codon]['f']
            target_f = usage_target.usage_table[aa][orig_codon]['f']

            df = abs(origin_f - target_f)
            new_codon = orig_codon

            codon['origin_f'] = origin_f
            codon['target_f'] = target_f
            codon['initial_df'] = df

            for item in usage_target.usage_table[aa]:
                if item != orig_codon:
                    new_df = abs(origin_f - usage_target.usage_table[aa][item]['f'])
                    if new_df < df:
                        df = new_df
                        new_codon = item

            codon['final_df'] = df
            codon['new'] = new_codon

    def construct_new_sequence(self):
        tmp = []
        for codon in self.codons:
            tmp.append(codon['new'])

        self.harmonized_sequence = Seq(''.join(tmp), IUPAC.unambiguous_dna)
        self.harmonized_translated_sequence = self.translate_sequence(self.harmonized_sequence, cds=True)

    def verify_harmonized_sequence(self):
        if str(self.original_translated_sequence) == str(self.harmonized_translated_sequence):
            return True
        else:
            return False

#    def align_sequences(self, aligner='tcoffee'):
#        if aligner == 'clustalw':
#        elif aligner == 'clustalo':
#        elif aligner == 'tcoffee':
#            from Bio.Align.Applications import TCoffeeCommandline

#            cline = TCoffeeCommandline()
#        elif aligner == 'mafft':
#        elif aligner == 'muscle':
#        else:
#            raise NameError('Invalid alignment tool specified: {}'.format(aligner))


def chunks(string, n):
    """Produce n-character chunks from string s."""
    for start in range(0, len(string), n):
        yield string[start:start + n]


#def main():
#    print("CHarm {}".format(__version__))
#    seq = Sequence("ATGTGCTAA")


#if __name__ == "__main__":
#    main()

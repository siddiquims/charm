from urllib.request import Request, urlopen
from urllib.error import URLError

try:
    # BeautifulSoup is used for parsing the HTML formatted codon usage tables on 'http://www.kazusa.or.jp/codon/'
    from bs4 import BeautifulSoup
except ImportError as e:
    print('ERROR: {}'.format(e.msg))
    exit(1)


class CodonUsageTable():
    """
    Provides a representation of a specific codon usage table
    url             - String; URL from which the usage table can be obtained.
                      Usually http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=<id>&aa=<num>&style=N
    use_frequency   - Boolean; Defines whether usage frequencies/1000 are used instead of fractions. Defaults to 'False'
    """

    def __init__(self, url, use_frequency=False):
        self.url = url
        self.usage_table = {}
        self.use_frequency = use_frequency
        self.fetch_codon_usage_table()

    def add_to_table(self, codon, aa, frequency):
        """
        Add codon and usage frequency to table
        codon     - String; e.g. 'ATG'
        aa        - String; Corresponding amino acid (e.g. 'M')
        frequency - Float: Usage fraction or frequency/1000
        """

        if aa in self.usage_table:
        # If the aa is already present in the table, just add the new codon
            self.usage_table[aa][codon] = {'f': frequency}
        else:
            # Else, create a new entry for the aa and add the codon
            self.usage_table[aa] = {}
            self.usage_table[aa][codon] = {'f': frequency}

    def fetch_codon_usage_table(self):
        """
        Fetch the codon table from http://www.kazusa.or.jp/codon
        """

        request = Request(self.url)
        try:
            # attempt to recieve the html file from the server
            opener = urlopen(request)
            response = opener.read()
            soup = BeautifulSoup(response)
        except URLError as e:
            # if this fails, print an error and exit
            if hasattr(e, 'reason'):
                print('Failed to reach server: %s' % e.reason)
            if hasattr(e, 'code'):
                print('Server responded with HTTP error code: %s' % e.code)
            exit(1)
            # otherwise read the file and hand over to BeautifulSoup for parsing

        # Parse the HTML response and look for a <pre></pre> section containing the usage table
        table_string = str(soup.pre)

        table_string = table_string.replace('<pre>\n', '')  # remove <pre></pre> tags
        table_string = table_string.replace('\n</pre>', '')

        table_lines = table_string.split('\n')  # Split in lines at the linebreak '\n'
        for line in table_lines:  # Iterate over the lines
            lines = line.split(')')  # Splitting the lines at ")" will result in substrings representing a
            # a single codon each
            for codon_raw in lines:
                codon_raw = codon_raw.strip()  # strip whitespace characters from the substring

                codon = codon_raw[:3].strip().replace('U', 'T')  # The first three characters are the codon
                if codon:
                    aa = codon_raw[4:5].strip()  # Position 5 is the aa in one letter code
                    fraction = float(codon_raw[6:10].strip())  # position 6 to 10 is the fraction;

                    frequency = float(codon_raw[11:15].strip())  # position 11 to 14 is the usage frequency/1000
                    # convert to float
                    if self.use_frequency:
                        # add either frequency or fraction to the codon table
                        self.add_to_table(codon, aa, frequency)
                    else:
                        self.add_to_table(codon, aa, fraction)
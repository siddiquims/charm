class CodonUsageTable():
    def __init__(self, url, use_frequency=False):
        self.url = url
        self.usage_table = {}
        self.use_frequency = use_frequency
        self.fetch_codon_usage_table()


    def add_to_table(self, codon, aa, frequency):
        """
        Add codon and usage frequency to table
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

        # Parse the HTML response and look for <pre></pre> section, containing the usage table
        table_string = str(soup.pre)

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
                    fraction = float(codon_raw[6:10].strip()) # position 6 to 10 is the fraction;

                    frequency = float(codon_raw[11:15].strip()) # position 11 to 14 is the usage frequency/1000
                    # convert to float
                    if self.use_frequency:
                        self.add_to_table(codon, aa, frequency)
                    else:
                        self.add_to_table(codon, aa, fraction)
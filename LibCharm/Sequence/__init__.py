from operator import itemgetter

try:
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    from Bio.Data import CodonTable
    from Bio.Data import IUPACData
except ImportError as e:
    print('ERROR: {}'.format(e.msg))
    exit(1)

from ..CodonUsageTable import CodonUsageTable


class Sequence():
    """
    Provides methods for storage and manipulation of sequences
    """

    def __init__(self, sequence, origin_id, host_id, translation_table_origin=1, translation_table_host=1,
                 use_frequency=False, lower_threshold=None, strong_stop=True, lower_alternative=True,
                 use_replacement_table=True, use_highest_frequency_if_ambiguous=True):
        """
        Initialize the Sequence object
        sequence                    - DNA or RNA sequence as Bio.Seq object or string. This can for example be
                                      generated by using BioPython directly or by loading a FASTA file using
                                      LibCharm.IO.load_file
        origin_id                   - Species id of the origin organism (can be found in the URL at
                                      http://www.kazusa.or.jp/codon)
        host_id                     - Species id of the host organism (can be found in the URL at
                                      http://www.kazusa.or.jp/codon)
        translation_table_host      - Integer; Genetic code used by the target host organism. Corrensponds to one of
                                      the translation tables listed here:
                                      http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        translation_table_origin    - Integer; Genetic code used by the target host organism. Corrensponds to one of
                                      the translation tables listed here:
                                      http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        use_frequency               - Boolean; Use frequency per thousand instead of fraction during the assessment of
                                      the codon usage
        lower_threshold             - Float; Threshold that defines the minimum codon usage that is considered
                                      appropriate. By default, a harmonized codon can only be lower than this
                                      threshold, if the original codon usage in the original organism is lower than
                                      this threshold, too.
        strong_stop                 - Boolean; Defines whether a strong stop codon (e.g. TAA in bacterial hosts) should
                                      be used. This may cause the stop codon not to be perfectly harmonized.
        lower_alternative           - Boolean; Defines whether the lower or higher usage codon should be used if df for
                                      two alternative codons is equal.
        use_replacement_table       - Boolean; If true, do not compute the harmonization for every single codon in the
                                      sequence, but for every unique codon in the sequence. This is done by default as
                                      it is much faster.
        use_highest_frequency_if_ambiguous - If the sequence contains ambigous codons (e.g. GCN),
        """

        # setting threshold if provided, otherwise fall back to defaults
        if not lower_threshold:
            if use_frequency:
                if not lower_threshold:
                    self.lower_threshold = 5
            else:
                if not lower_threshold:
                    self.lower_threshold = 0.1
        else:
            self.lower_threshold = lower_threshold

        # set other variables to provided values or defaults
        self.strong_stop = strong_stop
        self.lower_alternative = lower_alternative
        self.use_replacement_table = use_replacement_table
        self.use_frequency = use_frequency
        self.use_highest_frequency_if_ambiguous = use_highest_frequency_if_ambiguous

        self.ambiguous_dna_letters = list(set(IUPACData.ambiguous_dna_letters) - set(IUPACData.unambiguous_dna_letters))

        # check if translation table id is > 15. Values > 15 cannot be mapped to http://www.kazusa.or.jp/codon/!
        if translation_table_origin > 15 or translation_table_host > 15:
            raise ValueError('Though the NCBI lists more than 15 translation tables, CHarm is limited to the '
                             'first 15 as listed on \'http://www.kazusa.or.jp/codon/\'.')
        # Set translation table for original sequence
        self.translation_table_origin = CodonTable.ambiguous_dna_by_id[int(translation_table_origin)]
        self.translation_table_host = CodonTable.ambiguous_dna_by_id[int(translation_table_host)]
        # Reformat and sanitize sequence string (remove whitespaces, change to uppercase)
        if type(sequence) is 'str':
            try:
                # if a string is provided, check if it contains U and not T to distinguish between RNA and DNA
                if 'U' in sequence and not 'T' in sequence:
                    seq = Seq(''.join(sequence.upper().split()), IUPAC.ambiguous_rna)
                    # if RNA, convert to DNA alphabet
                    self.original_sequence = seq.back_transcribe()
                else:
                    self.original_sequence = Seq(''.join(sequence.upper().split()), IUPAC.ambiguous_dna)
            except ValueError as e:
                print('ERROR: {}'.format(e))
                exit(1)
        else:
            self.original_sequence = sequence
        self.original_translated_sequence = self.translate_sequence(self.original_sequence,
                                                                    self.translation_table_origin, cds=True)
        self.harmonized_sequence = ''

        self.usage_origin = CodonUsageTable('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?'
                                            'species={}&aa={}&style=N'.format(origin_id,
                                                                              translation_table_origin),
                                            self.use_frequency)
        self.usage_host = CodonUsageTable('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?'
                                          'species={}&aa={}&style=N'.format(host_id,
                                                                            translation_table_host),
                                          self.use_frequency)
        self.codons = self.split_original_sequence_to_codons()
        self.harmonize_codons()
        self.harmonized_sequence = self.construct_new_sequence()
        self.harmonized_translated_sequence = self.translate_sequence(self.harmonized_sequence,
                                                                      self.translation_table_host, cds=True)

    @staticmethod
    def chunks(string, n):
        """
        Produce n-character chunks from string.
        string  - string to be sliced
        n       - number of characters per chunk
        """
        for start in range(0, len(string), n):
            yield string[start:start + n]

    def translate_sequence(self, sequence, translation_table, cds=True, to_stop=False):
        """
        Translate a given DNA or RNA sequence into an amino acid sequence.

        sequence - The input sequence as Bio.Seq object
        cds      - Whether the input sequence is a coding region or not
        to_stop  - Only translate up to the first stop codon
        """

        translated_sequence = None
        try:
            translated_sequence = sequence.translate(table=translation_table, cds=cds, to_stop=to_stop)
        except CodonTable.TranslationError as e:
            print("Error during translation: ", e)
            print("This might be just fine if an additional stop codon was found at the end of the sequence.")
            return self.translate_sequence(sequence, translation_table, cds=False, to_stop=True)
        except KeyError as e:
            print("Error during translation: ", e)
            exit(1)
        return translated_sequence

    def get_harmonized_codons(self):
        """
        Returns a list of all harmonized codons
        """
        harmonized_codons = []
        for codon in self.codons:
            if str(codon['original']) != str(codon['new']) or codon['ambiguous']:
                harmonized_codons.append(codon)

        return harmonized_codons

    def split_original_sequence_to_codons(self):
        """
        Splits the sequence into codons.

        returns list 'codons' with
        position   - position of codon in the sequence (1 ... end)
        original   - original codon
        new        - new codon after harmonization
        origin_f   - usage frequency/fraction of the original codon in the origin organism
        target_f   - usage frequency/fraction of the new codon in the target host
        initial_df - difference in usage frequency/fraction of the original codon between origin organism and target
                     host
        final_df   - difference in usage frequency/fraction of the original codon in the origin organism and the
                     new codon after harmonization in the target host
        aa         - amino acid coded by the codon
        """

        codons = []
        position = 0

        for codon in self.chunks(self.original_sequence, 3):
            position += 1
            aa = str(codon.translate(table=self.translation_table_origin))

            codons.append({'position': int(position),
                           'original': str(codon),
                           'ambiguous': None,
                           'new': None,
                           'origin_f': None,
                           'target_f': None,
                           'initial_df': None,
                           'final_df': None,
                           'aa': aa})
        return codons

    def choose_wobble_codon(self, usage_table, codon, aa, highest_f):

        #codons = usage_table[aa].keys()
        codons = []
        for k in IUPACData.ambiguous_dna_values.keys():
            if k in self.ambiguous_dna_letters:
                values = list(IUPACData.ambiguous_dna_values[k])
                for v in values:
                    if k in codon:
                        codons.append(codon.replace(k, v))

        max_f = [None, 0.0]
        min_f = [None, 100.0]
        for c in codons:
            f = usage_table[aa][c]['f']
            if f > max_f[1]:
                max_f = [c, f]
            if f < min_f[1]:
                min_f = [c, f]

        if highest_f:
            return max_f
        else:
            return min_f

    def sort_replacement_codons(self, codons):

        for codon in codons:

            aa = codon['aa']
            orig_codon = str(codon['original'])

            if str(codon['original'])[-1] in self.ambiguous_dna_letters:
                orig_unambiguous_codon = self.choose_wobble_codon(self.usage_origin.usage_table, orig_codon, aa,
                                                                  self.use_highest_frequency_if_ambiguous)[0]
                codon['ambiguous'] = True
            else:
                orig_unambiguous_codon = orig_codon

            origin_f = self.usage_origin.usage_table[aa][orig_unambiguous_codon]['f']
            target_f = self.usage_host.usage_table[aa][orig_unambiguous_codon]['f']

            df = abs(origin_f - target_f)

            codon['origin_f'] = origin_f
            codon['target_f'] = target_f
            codon['initial_df'] = df

            codon_substitutions = []
            stop_codons = []

            for item in self.usage_host.usage_table[aa]:

                f_target_new = self.usage_host.usage_table[aa][item]['f']
                df_new = abs(origin_f - f_target_new)

                if aa == '*' and self.strong_stop:
                    stop_codons.append((item, df_new, f_target_new))
                else:
                    if f_target_new < self.lower_threshold < origin_f:
                        add = False
                    else:
                        if df_new < df:
                            add = True
                        else:  # target_f == 0:
                            add = True

                    if add:
                        codon_substitutions.append((item, df_new, f_target_new))

            if codon_substitutions:
                # sort the possible substitutions by df and frequency in target host
                sorted_codon_substitutions = sorted(codon_substitutions, key=itemgetter(1, 2))

                chosen_codon_index = 0  # choose lowest df by default

                if len(codon_substitutions) >= 2:

                    if (sorted_codon_substitutions[0][1] == sorted_codon_substitutions[1][1]) and not \
                            (sorted_codon_substitutions[0][2] == sorted_codon_substitutions[1][2]):
                        # if df of the first possible substitutions are identical
                        if not self.lower_alternative and len(sorted_codon_substitutions) > 1:
                            # choose the one with the higher frequency in target host if lower_alternative == False
                            chosen_codon_index += 1

                codon['final_df'] = sorted_codon_substitutions[chosen_codon_index][1]
                codon['target_f'] = sorted_codon_substitutions[chosen_codon_index][2]
                codon['new'] = sorted_codon_substitutions[chosen_codon_index][0]

            else:
                if aa == '*' and self.strong_stop:  # if this is a stop codon and we want a strong stop codon
                    sorted_stop_codons = sorted(stop_codons, key=itemgetter(2))  # sort by frequency in target host

                    # choose the codon with the highest usage frequency
                    codon['final_df'] = sorted_stop_codons[-1][1]
                    codon['target_f'] = sorted_stop_codons[-1][2]
                    codon['new'] = sorted_stop_codons[-1][0]
                else:
                    # if nothing fits better, leave the original codon in place
                    codon['final_df'] = df
                    codon['new'] = orig_unambiguous_codon
        return codons

    def compute_replacement_table(self):
        """
        Generates a list of unique codons and harmonize their codon usage. This list is returned and can be used
        to replace codons in a much longer list without the need to compute the codon substitution for every single
        position.
        """
        unique_codons = []
        unique_codons_triplets = []

        for codon in self.codons:
            if len(unique_codons) == 0 or codon['original'] not in unique_codons_triplets:
                unique_codons_triplets.append(codon['original'])
                unique_codons.append(codon)

        return self.sort_replacement_codons(unique_codons)

    def harmonize_codons(self):
        """
        Harmonizes the codon usage of self.original_sequence. This can either be done per codon or by
        computing a replacement table first (default). The second approach is much faster for long sequences but
        not as flexible.
        """

        if self.use_replacement_table:
            # This is a much faster approach, but not as flexible as the substitution is only done per codon and cannot
            # be expanded to its surroundings.
            codon_substitutions = self.compute_replacement_table()
            for codon in self.codons:
                for new_codon in codon_substitutions:
                    if codon['original'] == new_codon['original']:
                        for key, value in new_codon.items():
                            if key != 'position':
                                codon[key] = value

        else:
            self.codons = self.sort_replacement_codons(self.codons)

        return self.codons

    def construct_new_sequence(self):
        """
        Constructs the harmonized sequence out of the original and substituted codons
        in self.codons
        """
        tmp = []
        for codon in self.codons:
            tmp.append(codon['new'])

        harmonized_sequence = Seq(''.join(tmp), IUPAC.unambiguous_dna)
        #harmonized_translated_sequence = self.translate_sequence(harmonized_sequence, cds=True)

        return harmonized_sequence

    def verify_harmonized_sequence(self):
        """
        Verifies that the translation of the original and harmonized sequence is identical.
        This has to be true, but might fail due to potential errors in the algorithm.
        """
        if str(self.original_translated_sequence) == str(self.harmonized_translated_sequence):
            return True
        else:
            return False

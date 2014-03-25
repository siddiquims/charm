"""Handling data input from sequence files in multiple formats"""
try:
    from Bio import SeqIO
    from Bio.Data import CodonTable
    from Bio.Alphabet import IUPAC
except ImportError as e:
    print('ERROR: {}'.format(e.msg))
    exit(1)


def load_file(filename, file_format="fasta"):
    """
    Load sequence from file and returns sequence as Bio.Seq object
    :param filename:     String; Path and filename of input sequence file
    :param file_format:  String; Format to be used. Refer to Biopython docs for available formats. Defaults to 'fasta'
    """
    content = None
    try:
        # assume sequence is DNA
        content = SeqIO.read(filename, file_format, IUPAC.ambiguous_dna)
    except ValueError as error:
        # if this fails, try RNA instead
        print('ERROR: {}'.format(error))
        try:
            content = SeqIO.read(filename, file_format, IUPAC.ambiguous_rna)
        except ValueError as error:
            # if this fails, too, raise exception and exit with error code 1
            print('ERROR: {}'.format(error))
            exit(1)

    # if some kind of data could be read, return the sequence object
    if content:
        seq = content.seq
        return seq
    # else return None
    else:
        return None

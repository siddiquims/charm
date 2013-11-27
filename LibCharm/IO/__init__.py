try:
    from Bio import SeqIO
    from Bio.Alphabet import IUPAC
except ImportError as e:
    print('ERROR: {}'.format(e.msg))
    exit(1)

def load_file(filename, file_format="fasta"):
    """
    Load sequence from file in FASTA file_format
    filename - Path and filename of input sequence file
    """
    content = None
    try:
        content = SeqIO.read(filename, file_format, IUPAC.unambiguous_dna)
    except ValueError as e:
        print('ERROR: {}'.format(e.msg))
        try:
            content = SeqIO.read(filename, file_format, IUPAC.unambiguous_rna)
        except ValueError as e:
            print('ERROR: {}'.format(e.msg))
            exit(1)

    if content:
        seq = content.seq
        return seq
    else:
        return None

def file_open(filename, format="fasta"):
    """
    Load sequence from file in FASTA format
    filename - Path and filename of input sequence file
    """
    try:
        content = SeqIO.read(filename, format, IUPAC.unambiguous_dna)
    except:
        try:
            content = SeqIO.read(filename, format, IUPAC.unambiguous_rna)
        except:
            exit(1)

    seq = content.seq
    return seq

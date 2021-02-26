import pandas


def align_pair(
        query_seq,
        reference_seq,
        local=False,
        gap_open_penalty=11,
        gap_extension_penality=1,
        substitution_matrix="blosum62",
        alignment_function=None):
    """
    Align two protein or DNA sequences.

    By default, a protein substitution matrix (blosum62) is used. If you are
    aligning DNA or RNA, you should use a nucleotide substitution matrix by
    passing, for example, substitution_matrix="dnafull".

    This is a thin wrapper over the Parasail library implementation.

    Returns a pandas.Series with the results of the alignment.

    Parameters
    ----------
    query_seq : string
        First sequence to align
    reference_seq : string
        Second sequence to align.
    local : boolean
        If True, a local alignment is performed using the Smith-Waterman
        algorithm. This means that gaps at the beginning or end of the
        sequences are not penalized, and only the part of the sequences that
        align are returned.

        If False, a global alignment is performed using the Needleman-Wunsch
        algorithm. This means that the two sequences will be aligned in their
        entirety.
    gap_open_penalty : int
        Penality for starting a gap
    gap_extension_penality : int
        Penalty for extending a gap
    substitution_matrix : string
        Name of substitution matrix. Examples: "blosum62", "blosum90",
        "dnafull", "pam100". If you are aligning DNA or RNA you should use a
        nucleotide substitution matrix, such as "dnafull".

        Full list of supported matrices:
        https://github.com/jeffdaily/parasail/tree/master/parasail/matrices
    alignment_function : function
        Advanced use. If you know the underlying parasail alignment function
        you would like to use, you can pass it here. Otherwise a reasonable
        default is used.

    Returns
    -------
    pandas.Series with keys:
        query : string
            Aligned query sequence
        reference : string
            Aligned reference sequence
        correspondence : string
            Characters (similar to BLAST "midline") indicating the
            correspondence between query and reference strings.
        score : int
            Alignment score. Higher indicates a better alignment.
    """
    import parasail

    matrix = parasail.Matrix(substitution_matrix)

    if alignment_function is None:
        if local:
            alignment_function = parasail.sw_trace_striped_16
        else:
            alignment_function = parasail.nw_trace_striped_16

    raw_results = alignment_function(
        query_seq,
        reference_seq,
        gap_open_penalty,
        gap_extension_penality,
        matrix)

    tb = raw_results.get_traceback()
    result = pandas.Series(dtype=object)
    result["query"] = tb.query
    result["reference"] = tb.ref
    result["correspondence"] = tb.comp
    result["score"] = raw_results.score
    return result
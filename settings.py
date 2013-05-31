default_sample_settings = {
    'minimum_length': 20,       # minimum tag length after adapter trimming
    'quality_threshold': 25,    # okay if `quality_percentage'% of bases in a read are
    'quality_percentage': 90,   #   as good as `quality_threshold' or better.
    'piranha_bin_size': 30,

    # Default parameters for clustered error-based binding site identification.
    # These parameters are used only for showing motifs and secondary structure on the
    # default result report. The other data files are free from these settings
    # as they cover all possible methods and broad range of significance levels.
    'default_error_scoring': 'entropy',
    'default_error_score_cutoff': 0.8,
    'default_error_depth': 50,
}

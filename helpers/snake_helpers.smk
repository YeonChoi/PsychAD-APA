
def read_table(path):
    from pathlib import Path
    READERS = {
        '.parquet': pd.read_parquet,
        '.csv':     pd.read_csv,
        '.csv.gz':     pd.read_csv,
        '.tsv':     lambda p: pd.read_csv(p, sep='\t'),
        '.tsv.gz':     lambda p: pd.read_csv(p, sep='\t'),
        '.txt':     lambda p: pd.read_csv(p, sep='\t'),
        '.txt.gz':     lambda p: pd.read_csv(p, sep='\t'),
        '.feather': pd.read_feather,
        '.pkl':     pd.read_pickle,
    }

    path = Path(path)
    suffixes = [s.lower() for s in path.suffixes]

    ext = suffixes[-2] if suffixes and suffixes[-1] in {'.gz', '.bz2', '.zst', '.xz'} else (suffixes[-1] if suffixes else '')

    if ext not in READERS:
        raise ValueError(f'Unsupported file extension: {ext} ({path})')
    if not path.exists():
        raise FileNotFoundError(path)

    return READERS[ext](path)

import pandas as pd


def find_chr_col(df: pd.DataFrame) -> str:
    if "Chromosome" in df.columns:
        return "Chromosome"
    else:
        return "chr"


def find_start_col(df: pd.DataFrame) -> str:
    if "Start" in df.columns:
        return "Start"
    else:
        return "start"


def find_end_col(df: pd.DataFrame) -> str:
    if "End" in df.columns:
        return "End"
    else:
        return "end"


def find_strand_col(df: pd.DataFrame) -> str:
    if "Strand" in df.columns:
        return "Strand"
    else:
        return "strand"

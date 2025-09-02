import pandas as pd

import sqlite3
from .genes import SEP
from .genomic.location import NA

BLACKLIST_QUERY = f"""
SELECT DISTINCT 
    chr, start, end notes 
FROM blacklist
WHERE 
    chr = :chromosome AND 
    :start <= end AND :end >= start
ORDER BY chr, start;
"""


class BlacklistAnnotation:
    def __init__(self, blacklist_db_file: str):
        self._blacklist_db_file = blacklist_db_file

    def annotate_blacklist(self, df_query: pd.DataFrame):
        conn = sqlite3.connect(self._blacklist_db_file)

        print(self._blacklist_db_file)

        # Create a cursor object
        cursor = conn.cursor()

        blacklists = []

        for _, row in df_query.iterrows():
            records = []

            cursor.execute(
                BLACKLIST_QUERY,
                {
                    "chromosome": row["Chromosome"],
                    "start": row["Start"],
                    "end": row["End"],
                },
            )

            # Fetch all records
            records.extend(cursor.fetchall())

            if len(records) > 0:
                blacklists.append(
                    SEP.join([f"{rec[0]}:{rec[1]}-{rec[2]}" for rec in records])
                )
            else:
                blacklists.append(NA)

        cursor.execute("SELECT genome, version FROM info")
        genome_info = cursor.fetchone()
        version = genome_info[1]

        df_query[version] = blacklists

        cursor.close()
        conn.close()

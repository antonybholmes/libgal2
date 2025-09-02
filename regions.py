import collections
import pandas as pd

import sqlite3


from .genes import SEP
from .genomic import NA
from . import utils, genomic

# BLACKLIST_QUERY = f"""
# SELECT DISTINCT
#     chr, start, end, notes
# FROM regions
# WHERE
#     chr = :chromosome AND
#     :start <= end AND :end >= start
# ORDER BY chr, start;
# """

REGION_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    r.start,
    r.end
FROM query_regions q
JOIN regions r ON 
    r.chr = q.chr AND
    r.start <= q.end AND 
    r.end >= q.start
ORDER BY q.location;
"""

TEMP_QUERY_TABLE_SQL = f"""
    CREATE TEMP TABLE IF NOT EXISTS query_regions (
    location TEXT NOT NULL,
    chr TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand TEXT NOT NULL DEFAULT '+'
);
"""

INSERT_TEMP_QUERY = f"""
    INSERT INTO query_regions (location, chr, start, end)
    VALUES (:location, :chr, :start, :end)
"""

TEMP_INDEX_QUERY_TABLE_REGION_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, start, end)"
)

DELETE_QUERY_TABLE_SQL = f"""
    DELETE FROM query_regions
"""


def row_to_dict(row):
    location = row[0]
    chr = row[1]
    start = row[2]
    end = row[3]

    annotation = {
        "location": location,
        "chr": chr,
        "start": start,
        "end": end,
    }

    return annotation


class RegionAnnotation:
    def __init__(self):
        self._tad_db_file = None
        self._conn = None
        self._cursor = None

    def open(
        self,
        db: str,
    ):
        #  close any existing connections
        self.close()

        self._db = db

        print(f"Opening database connection to {self._db}")

        self._conn = sqlite3.connect(self._db)

        # Create a cursor object
        self._cursor = self._conn.cursor()

        self._cursor.execute(TEMP_QUERY_TABLE_SQL)
        self._cursor.execute(TEMP_INDEX_QUERY_TABLE_REGION_SQL)

    def close(self):
        if self._cursor:
            self._cursor.close()
        if self._conn:
            self._conn.close()

    def annotate_df(self, df_query: pd.DataFrame, header: str = "Regions"):
        print(f"Adding {header} annotations to dataframe...")
        locs = []

        chr_col = utils.find_chr_col(df_query)
        start_col = utils.find_start_col(df_query)
        end_col = utils.find_end_col(df_query)

        for _, row in df_query.iterrows():
            loc = genomic.Location(
                chr=row[chr_col],
                start=row[start_col],
                end=row[end_col],
                strand="+",
            )
            locs.append(loc)

        annotations = self.annotate(locs)

        print(annotations)

        df_query[header] = [
            (SEP.join([str(l) for l in a]) if len(a) > 0 else genomic.NA)
            for a in annotations
        ]

        print("Done.")

    def annotate(self, locations: list[genomic.Location]):
        self._cursor.execute(DELETE_QUERY_TABLE_SQL)

        print(f"Annotating regions using {self._db}")

        queries = []
        for loc in locations:
            queries.append(
                {
                    "location": str(loc),
                    "chr": loc.chr,
                    "start": loc.start,
                    "end": loc.end,
                    "strand": loc.strand,
                }
            )

        self._cursor.executemany(
            INSERT_TEMP_QUERY,
            queries,
        )

        self._cursor.execute(REGION_QUERY)

        annotation_map = collections.defaultdict(set)

        for c in self._cursor:
            d = row_to_dict(c)

            chr = d["chr"]
            start = d["start"]
            end = d["end"]
            loc = genomic.location.parse_location(d["location"])
            region_loc = genomic.Location(chr, start, end)

            annotation_map[loc].add(region_loc)

        ret = []

        for loc in locations:
            ret.append(list(sorted(annotation_map[loc])))

        print("Done.")

        return ret


class BlacklistAnnotation(RegionAnnotation):
    def open(self, blacklist_db_file: str):
        super().open(blacklist_db_file)

    def annotate_df(self, df_query: pd.DataFrame, header: str = "Blacklist regions"):
        super().annotate_df(df_query, header)


class CentromereAnnotation(RegionAnnotation):
    def open(self, centromere_db_file: str):
        super().open(centromere_db_file)

    def annotate_df(self, df_query: pd.DataFrame, header: str = "Centromere regions"):
        super().annotate_df(df_query, header)


class TelomereAnnotation(RegionAnnotation):
    def open(self, telomere_db_file: str):
        super().open(telomere_db_file)

    def annotate_df(self, df_query: pd.DataFrame, header: str = "Telomere regions"):
        super().annotate_df(df_query, header)

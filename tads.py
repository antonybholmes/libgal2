import collections
import json
import pandas as pd

import sqlite3

from . import genomic, utils
from .genes import SEP


TAD_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    t.start,
    t.end,
    t.gene_ids, 
    t.gene_names
FROM query_regions q
JOIN tads t ON 
    t.chr = q.chr AND
    t.start <= q.end AND 
    t.end >= q.start
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
    INSERT INTO query_regions (location, chr, start, end, strand)
    VALUES (:location, :chr, :start, :end, :strand)
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
    gene_ids = row[4]
    gene_names = row[5]

    annotation = {
        "location": location,
        "chr": chr,
        "start": start,
        "end": end,
        "gene_ids": gene_ids,
        "gene_names": gene_names,
    }

    return annotation


class TADAnnotation:
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

    def annotate_df(self, df_query: pd.DataFrame):
        print("Adding TAD annotations to dataframe")
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

        if len(annotations) > 0:
            df_query["TAD domains"] = [
                SEP.join([str(l["location"]) for l in a["tads"]]) for a in annotations
            ]

            df_query["Genes in same TAD domain"] = [
                SEP.join(
                    [
                        ",".join([str(g["gene_name"]) for g in l["annotations"]])
                        for l in a["tads"]
                    ]
                )
                for a in annotations
            ]
        else:
            df_query["TAD domains"] = genomic.NA
            df_query["Genes in same TAD domain"] = genomic.NA

    def annotate(self, locations: list[genomic.Location]) -> list[dict]:
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

        self._cursor.execute(TAD_QUERY)

        annotation_map = collections.defaultdict(lambda: collections.defaultdict(set))

        for c in self._cursor:
            d = row_to_dict(c)

            chr = d["chr"]
            start = d["start"]
            end = d["end"]
            loc = genomic.location.parse_location(d["location"])
            tad_loc = genomic.Location(chr, start, end)
            ids = d["gene_ids"].split(",")
            names = d["gene_names"].split(",")

            ann = [{"gene_id": i, "gene_name": n} for i, n in zip(ids, names)]

            for a in ann:
                annotation_map[loc][tad_loc].add(json.dumps(a))

        ret = []

        for loc in locations:
            loc_entry = {"location": loc, "tads": []}

            for tad_loc in sorted(annotation_map[loc]):
                ann_set = annotation_map[loc][tad_loc]

                loc_entry["tads"].append(
                    {
                        "location": tad_loc,
                        "annotations": sorted(
                            [json.loads(d) for d in ann_set],
                            key=lambda x: x["gene_name"],
                        ),
                    }
                )

            ret.append(loc_entry)

        return ret

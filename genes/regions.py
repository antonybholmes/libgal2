import collections
from . import dataframe  # local module
import sqlite3
from .. import genomic
import json


GENE_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    g.start,
    g.end,
    g.gene_id, 
    g.gene_name,
    g.transcript_id,
    g.strand,
    g.feature
FROM query_regions q
JOIN gtf g ON 
    g.feature = :feature AND 
    g.seqname = q.chr AND
    g.start <= q.end AND 
    g.end >= q.start
ORDER BY q.location;
"""


def row_to_dict(row):
    location = row[0]
    chr = row[1]
    start = row[2]
    end = row[3]
    gene_id = row[4]
    gene_name = row[5]
    transcript_id = row[6]
    gene_strand = row[7]
    feature = row[8]

    annotation = {
        "location": location,
        "chr": chr,
        "start": start,
        "end": end,
        "transcript_id": transcript_id,
        "gene_id": gene_id,
        "gene_name": gene_name,
        "feature": feature,
        "strand": gene_strand,
    }

    return annotation


class Annotation:
    def __init__(
        self,
        closest_n: int = 5,
        max_distance: int = 2000000,
        promoter_lim: list[int] = [2000, 1000],
    ):
        self._closest_n = closest_n
        self._max_distance = max_distance
        self._promoter_lim = promoter_lim
        self._conn = None
        self._cursor = None
        self._queries = []
        self._db = None

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

        self._cursor.execute(dataframe.TEMP_QUERY_TABLE_SQL)
        self._cursor.execute(dataframe.TEMP_INDEX_QUERY_TABLE_REGION_SQL)
        self._cursor.execute(dataframe.TEMP_INDEX_MID_SQL)

    def close(self):
        if self._cursor:
            self._cursor.close()
        if self._conn:
            self._conn.close()

    def annotate(self, locations: list[genomic.Location], feature: str = "gene"):
        self._cursor.execute(dataframe.DELETE_QUERY_TABLE_SQL)

        print(f"Annotating regions using {self._db}")

        queries = []
        for loc in locations:
            queries.append(
                {
                    "location": str(loc),
                    "chr": loc.chr,
                    "start": loc.start,
                    "end": loc.end,
                    "midpoint": loc.mid,
                    "strand": loc.strand,
                }
            )

        self._cursor.executemany(
            dataframe.INSERT_TEMP_QUERY,
            queries,
        )

        print(f"Processing {feature}s...")

        # find out which intronic regions are exonic

        self._cursor.execute(GENE_QUERY, {"feature": feature})

        annotation_map = collections.defaultdict(set)

        for c in self._cursor:
            d = row_to_dict(c)

            annotation_map[d["location"]].add(
                json.dumps(
                    {
                        "gene_id": d["gene_id"],
                        "gene_name": d["gene_name"],
                        "strand": d["strand"],
                    }
                )
            )

        return [
            {
                "location": loc,
                "annotations": sorted(
                    [json.loads(d) for d in annotation_map.get(str(loc), set())],
                    key=lambda x: x["gene_name"],
                ),
            }
            for loc in locations
        ]

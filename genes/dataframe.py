import collections
import json
import sys
from typing import Union
import pandas as pd
import pyranges as pr

import sqlite3


SEP = " | "
NA = "n/a"
PROMOTER_LIM = [2000, 1000]
TOP = 40

# Example: Execute a SELECT query
# TRANSCRIPT_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.gene_symbol,
#     transcript.transcript_id,
#     gene.chr,
#     transcript.tss,
#     gene.strand,
#     'intronic'
#     FROM gene
#     JOIN transcript ON gene.id = transcript.gene_id
#     WHERE gene.chr = :chromosome AND
#     :midpoint >= transcript.start AND :midpoint <= transcript.end
#     """

INTRONIC_JOIN_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    q.midpoint,
    g.gene_id, 
    g.gene_name,
    g.transcript_id,
    g.tss,
    g.strand,
    'intronic' AS type,
    ABS(q.midpoint - g.tss) as ab_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND 
    q.midpoint >= g.start AND 
    q.midpoint <= g.end
ORDER BY q.location;
"""

# EXON_QUERY = f"""SELECT DISTINCT g.gene_id,
#     g.gene_symbol,
#     transcript.transcript_id,
#     exon.tss,
#     gene.strand,
#     'exonic'
#     FROM gene, transcript, exon
#     WHERE gene.id = transcript.gene_id AND
#     transcript.id = exon.transcript_id AND
#     gene.chr = :chromosome AND
#     :midpoint >= exon.start AND :midpoint <= exon.end
#     """

EXON_JOIN_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    q.midpoint,
    g.gene_id, 
    g.gene_name,
    g.transcript_id,
    g.tss,
    g.strand,
    'exonic' AS type,
    ABS(q.midpoint - g.tss) as ab_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'exon' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint
ORDER BY q.location;
"""

IS_INTRONIC_QUERY = f"""
SELECT DISTINCT q.row_idx
FROM intronic_query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint AND
    g.transcript_id = q.transcript_id
ORDER BY q.row_idx;
"""

IS_EXONIC_QUERY = f"""
SELECT DISTINCT q.row_idx
FROM intronic_query_regions q
JOIN gtf g ON 
    g.feature = 'exon' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint AND
    g.transcript_id = q.transcript_id
ORDER BY q.row_idx;
"""

IS_PROMOTER_QUERY = f"""
SELECT DISTINCT q.row_idx
FROM intronic_query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND
    (g.tss - :promoter_lim_1) <= q.midpoint AND 
    (g.tss + :promoter_lim_2) >= q.midpoint AND 
    g.strand = :strand AND
    g.transcript_id = q.transcript_id
ORDER BY q.row_idx;
"""


# PROMOTER_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.gene_symbol,
#     transcript.transcript_id,
#     transcript.tss,
#     gene.strand,
#     'promoter'
#     FROM gene, transcript
#     WHERE
#     gene.strand = :strand AND
#     gene.chr = :chromosome AND
#     gene.id = transcript.gene_id AND
#     :midpoint >= (transcript.tss - :promoter_lim_1) AND :midpoint <= (transcript.tss + :promoter_lim_2)
#     """

PROMOTER_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    q.midpoint,
    g.gene_id, 
    g.gene_name,
    g.transcript_id, 
    g.tss,
    g.strand,
    'promoter' AS type,
    ABS(q.midpoint - g.tss) as ab_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND
    (g.tss - :promoter_lim_1) <= q.midpoint AND 
    (g.tss + :promoter_lim_2) >= q.midpoint AND 
    g.strand = :strand
ORDER BY q.location;
"""


# NEAREST_INTRONIC_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.gene_symbol,
#     transcript.transcript_id,
#     transcript.tss,
#     gene.strand,
#     ABS(transcript.tss - :midpoint) AS ab_dist,
#     'intronic'
#     FROM gene, transcript
#     WHERE
#     gene.chr = :chromosome AND
#     gene.id = transcript.gene_id
#     ORDER BY ab_dist, gene.gene_name
#     LIMIT 50
#     """

NEAREST_GENE_JOIN_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    q.midpoint,
    g.gene_id, 
    g.gene_name,
    g.transcript_id,
    g.tss,
    g.strand,
    'intergenic' AS type,
    g.start,
    g.end,
    ABS(q.midpoint - g.tss) as ab_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND
    g.seqname = q.chr
WHERE ab_dist < :max_distance
ORDER BY q.location, ab_dist, g.gene_name;
"""


# NEAREST_INTRONIC_JOIN_QUERY = f"""
# SELECT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.gene_name,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'intronic' AS type,
#     ABS(q.midpoint - g.tss) as ab_dist
# FROM query_regions q
# JOIN gtf g ON
#     g.seqname = q.chr AND
#     g.feature = 'transcript' AND
#     g.start <= q.midpoint AND
#     g.end >= q.midpoint
# WHERE ab_dist < :max_distance
# ORDER BY q.location, ab_dist, g.gene_name;
# """

NEAREST_EXON_JOIN_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    q.midpoint,
    g.gene_id, 
    g.gene_name,
    g.transcript_id, 
    g.tss,
    g.strand,
    'exonic' AS type,
    ABS(q.midpoint - g.tss) as ab_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'exon' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint
ORDER BY q.location, ab_dist, g.gene_name;
"""

NEAREST_PROMOTER_JOIN_QUERY = f"""
SELECT DISTINCT
    q.location,
    q.chr,
    q.midpoint,
    g.gene_id, 
    g.gene_name,
    g.transcript_id,
    g.tss,
    g.strand,
    'promoter' AS type,
    g.start,
    g.end,
    ABS(q.midpoint - g.tss) as ab_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND 
    (g.tss - :promoter_lim_1) <= q.midpoint AND 
    (g.tss + :promoter_lim_2) >= q.midpoint AND 
    g.strand = :strand
ORDER BY q.location, g.gene_name;
"""


TEMP_QUERY_TABLE_SQL = f"""
    CREATE TEMP TABLE IF NOT EXISTS query_regions (
    location TEXT NOT NULL,
    chr TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    midpoint INTEGER NOT NULL,
    strand TEXT NOT NULL DEFAULT '+'
);
"""

TEMP_INTRONIC_TABLE_SQL = f"""
    CREATE TEMP TABLE IF NOT EXISTS intronic_query_regions (
    row_idx INTEGER NOT NULL,
    location TEXT NOT NULL,
    chr TEXT NOT NULL,
    midpoint INTEGER NOT NULL,
    transcript_id TEXT NOT NULL
);
"""

INSERT_INTRONIC_QUERY = f"""
    INSERT INTO intronic_query_regions (row_idx, location, chr, midpoint, transcript_id)
    VALUES (:row_idx, :location, :chr, :midpoint, :transcript_id)
"""

DROP_INTRONIC_TABLE_SQL = f"""
    DROP TABLE IF EXISTS intronic_query_regions
"""

DELETE_INTRONIC_TABLE_SQL = f"""
    DELETE FROM intronic_query_regions
"""

DELETE_QUERY_TABLE_SQL = f"""
    DELETE FROM query_regions
"""

TEMP_INDEX_QUERY_TABLE_REGION_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, start, end)"
)

TEMP_INDEX_MID_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, midpoint)"
)

# TEMP_INTRONIC_TABLE_SQL = (
#     "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, start, end)"
# )

# TEMP_LOCATION_INDEX_TABLE_SQL = (
#     "CREATE INDEX IF NOT EXISTS idx_query_regions_location ON query_regions (location)"
# )

TEMP_INTRONIC_INDEX_SQL = "CREATE INDEX IF NOT EXISTS idx_intronic_query_regions ON intronic_query_regions (chr, midpoint)"

INSERT_TEMP_QUERY = f"""
    INSERT INTO query_regions (location, chr, start, end, midpoint, strand)
    VALUES (:location, :chr, :start, :end, :midpoint, :strand)
"""


def format_col_str(key: str, entrez_sorted_annotations: dict[str]) -> str:
    return SEP.join(str(ann[key]) for ann in entrez_sorted_annotations)


def format_col_set(
    key: str, entrez_sorted_annotations: dict[list[Union[str, int]]]
) -> str:
    return SEP.join(
        ",".join([str(x) for x in sorted(ann[key])])
        for ann in entrez_sorted_annotations
    )


def add_annotation_for_location_to_cols(
    annotations: list[dict], annotation_cols: list[list[Union[str, int]]]
):
    # check if contains intronic
    # is_intronic = len(list(filter(lambda x: x["type"] == "intronic", annotations))) > 0

    # if is_intronic:

    # if "intronic" in statuses:
    #    statuses.add("intragenic")

    # is_exonic = len(list(filter(lambda x: "exonic" in x["type"], annotations))) > 0

    # if is_exonic:
    # remove intronic annotations if also exonic as we prefer that
    #    annotations = list(filter(lambda x: "intronic" not in x["type"], annotations))

    # group by entrez
    gene_id_map = collections.defaultdict(dict)
    for ann in annotations:
        gene_id = ann["gene_id"]

        if gene_id not in gene_id_map:
            gene_id_map[gene_id] = {
                "gene_id": gene_id,
                "gene_symbol": "",
                "type": set(),
                "tss_distance": {"d": sys.maxsize, "transcript_id": ""},
                "strand": "+",
            }

        # entrez_map[ann["gene_id"]]["gene_id"].add(ann["gene_id"])

        gene_id_map[ann["gene_id"]]["gene_symbol"] = ann["gene_symbol"]
        # entrez_map[ann["gene_id"]]["transcript_id"] = ann["transcript_id"]
        gene_id_map[ann["gene_id"]]["type"].update(ann["type"])  # (ann["type"])
        # entrez_map[ann["gene_id"]]["tss_distance"] = ann["tss_distance"]

        # we use all annotations for labelling a region, but to reduce
        # noise, we keep the closest refeq by TSS dist per gene and report
        # only that so that we don't have long lists of refseqs for each
        # gene, which are likely redundant
        if abs(ann["tss_distance"]) < abs(
            gene_id_map[ann["gene_id"]]["tss_distance"]["d"]
        ):
            gene_id_map[ann["gene_id"]]["tss_distance"]["d"] = ann["tss_distance"]
            gene_id_map[ann["gene_id"]]["tss_distance"]["transcript_id"] = ann[
                "transcript_id"
            ]

        gene_id_map[ann["gene_id"]]["strand"] = ann["strand"]

    if len(gene_id_map) > 0:
        gene_id_sorted_annotations = []
        for gene_id in sorted(gene_id_map):
            # collapse tss distance
            ann = gene_id_map[gene_id]

            if "exonic" in ann["type"]:
                ann["type"].discard("intronic")

            # if ann["type"] contains intronic or exonic
            if "intronic" in ann["type"] or "exonic" in ann["type"]:
                ann["type"].add("intragenic")
                ann["type"].discard("intergenic")
            else:
                # ensure that it is marked as not being in a gene
                ann["type"].add("intergenic")

            gene_id_sorted_annotations.append(gene_id_map[gene_id])

        # transcript id
        annotation_cols[0].append(
            SEP.join(
                ann["tss_distance"]["transcript_id"]
                for ann in gene_id_sorted_annotations
            )
        )

        # gene id
        annotation_cols[1].append(format_col_str("gene_id", gene_id_sorted_annotations))

        annotation_cols[2].append(
            format_col_str("gene_symbol", gene_id_sorted_annotations)
        )

        annotation_cols[3].append(format_col_str("strand", gene_id_sorted_annotations))
        annotation_cols[4].append(
            SEP.join(
                str(ann["tss_distance"]["d"]) for ann in gene_id_sorted_annotations
            )
        )

        annotation_cols[5].append(format_col_set("type", gene_id_sorted_annotations))
    else:
        annotation_cols[0].append(NA)
        annotation_cols[1].append(NA)
        annotation_cols[2].append(NA)
        annotation_cols[3].append(NA)
        annotation_cols[4].append(NA)
        annotation_cols[5].append("intergenic")


def row_to_dict(row):
    location = row[0]
    chr = row[1]
    midpoint = row[2]
    gene_id = row[3]
    gene_symbol = row[4]
    transcript_id = row[5]
    transcript_tss = row[6]
    gene_strand = row[7]
    transcript_type = row[8]

    if gene_strand == "+":
        dist = midpoint - transcript_tss
    else:
        dist = transcript_tss - midpoint

    if dist == -145239:
        print(
            "Found matching distance:",
            location,
            dist,
            midpoint,
            transcript_tss,
            transcript_type,
        )

    annotation = {
        "location": location,
        "chr": chr,
        "midpoint": midpoint,
        "transcript_id": transcript_id,
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "type": [transcript_type],
        "tss_distance": dist,
        "strand": gene_strand,
    }

    return annotation


def row_to_annotation(row, annotation_map):
    # annotation = row_to_dict(row)

    annotation_map[row["location"]].add(
        json.dumps(row, sort_keys=True, separators=(",", ":"))
    )

    return row


# def row_to_dict_closest(row):
#     location = row[0]
#     chr = row[1]
#     midpoint = row[2]
#     gene_id = row[3]
#     gene_symbol = row[4]
#     transcript_id = row[5]
#     transcript_tss = row[6]
#     # ab_dist = row[7]
#     gene_strand = row[7]
#     annot_type = row[8]

#     dist = midpoint - transcript_tss

#     if gene_strand == "-":
#         dist = -dist

#     annotation = {
#         "location": location,
#         "chr": chr,
#         "midpoint": midpoint,
#         "transcript_id": transcript_id,
#         "gene_id": gene_id,
#         "gene_symbol": gene_symbol,
#         "type": annot_type,
#         "tss_distance": dist,
#         "strand": gene_strand,
#     }

#     return annotation


def row_to_closest_annotation(
    annotation: dict,
    closest_annotation_map: dict[dict[str]],
    used_symbols: dict[dict[int]],
):

    location = annotation["location"]
    gene_symbol = annotation["gene_symbol"]

    # if gene_symbol not in used_symbols[location]:
    #     if len(used_symbols[location]) < closest_n:
    #         used_symbols[location][gene_symbol] = len(used_symbols[location]) + 1

    # we keep the first 5 genes we encounter
    closest = used_symbols[location].get(gene_symbol, -1)

    # if closest == -1:
    # stop
    #    print("stopping")
    #    return True

    # print(gene_symbol, closest, closest_n)

    if closest != -1:
        closest_annotation_map[closest][location].add(
            json.dumps(annotation, sort_keys=True, separators=(",", ":"))
        )


class DataframeAnnotation:
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
        self._df_query = None
        self._prom_header = f"Relative To Gene (prom=-{promoter_lim[0]/1000}/+{promoter_lim[1]/1000} kb)"
        self._db = None

    def open(
        self,
        db: str,
        df_query: pd.DataFrame,
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
        self._cursor.execute(TEMP_INDEX_MID_SQL)

        self._queries = []
        for _, row in df_query.iterrows():
            location = f"{row['Chromosome']}:{row['Start']}-{row['End']}"
            midpoint = int((row["Start"] + row["End"]) / 2)
            self._queries.append(
                {
                    "location": location,
                    "chr": row["Chromosome"],
                    "start": row["Start"],
                    "end": row["End"],
                    "midpoint": midpoint,
                    "strand": "+",
                }
            )

        self._cursor.executemany(
            INSERT_TEMP_QUERY,
            self._queries,
        )

        self._cursor.execute(DROP_INTRONIC_TABLE_SQL)
        self._cursor.execute(TEMP_INTRONIC_TABLE_SQL)
        self._cursor.execute(TEMP_INTRONIC_INDEX_SQL)

        self._df_query = df_query

    def close(self):
        if self._cursor:
            self._cursor.close()
        if self._conn:
            self._conn.close()

    def annotate_genes(
        self,
    ):
        print(f"Annotating {len(self._df_query)} regions using {self._db}")

        annotation_map = collections.defaultdict(set)

        print("Processing introns...")

        self._cursor.execute(INTRONIC_JOIN_QUERY)

        # for row in cursor:
        #    row_to_annotation(row, annotation_map)

        rows = []
        for row in self._cursor:
            d = row_to_dict(row)

            rows.append(d)

        queries = []
        for idx, row in enumerate(rows):
            queries.append(
                {
                    "row_idx": idx,
                    "location": row["location"],
                    "chr": row["chr"],
                    "midpoint": row["midpoint"],
                    "transcript_id": row["transcript_id"],
                }
            )

        self._cursor.executemany(
            INSERT_INTRONIC_QUERY,
            queries,
        )

        print("Processing exons...")

        # find out which intronic regions are exonic

        self._cursor.execute(IS_EXONIC_QUERY)

        for c in self._cursor:
            if "exonic" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("exonic")

        for row in rows:
            # print(row)
            row_to_annotation(row, annotation_map)

        print("Processing promoters...")

        self._cursor.execute(
            PROMOTER_QUERY,
            {
                "strand": "+",
                "promoter_lim_1": self._promoter_lim[0],
                "promoter_lim_2": self._promoter_lim[1],
            },
        )

        for row in self._cursor:
            row_to_annotation(row_to_dict(row), annotation_map)

        self._cursor.execute(
            PROMOTER_QUERY,
            {
                "strand": "-",
                "promoter_lim_1": self._promoter_lim[1],
                "promoter_lim_2": self._promoter_lim[0],
            },
        )

        for row in self._cursor:
            row_to_annotation(row_to_dict(row), annotation_map)

        print("Adding annotations...")

        for row in self._cursor:
            row_to_annotation(row_to_dict(row), annotation_map)

        # if closest_n > 0:
        #     self.annotate_closest_genes(
        #         cursor, closest_n, max_distance, promoter_lim, closest_annotation_map
        #     )

        # add columns we are going to fill

        # self._df_query["Transcript Id"] = ""
        # self._df_query["Gene Id"] = ""
        # self._df_query["Gene Symbol"] = ""
        # self._df_query["Strand"] = ""
        # self._df_query["TSS Distance"] = ""
        # self._df_query[prom_header] = ""

        # if closest_n > 0:
        #     for i in range(1, closest_n + 1):
        #         df_query[f"#{i} Transcript Id"] = ""
        #         df_query[f"#{i} Gene Id"] = ""
        #         df_query[f"#{i} Gene Symbol"] = ""
        #         df_query[f"#{i} Strand"] = ""
        #         df_query[f"#{i} TSS Distance"] = ""
        #         df_query[f"#{i} {prom_header}"] = ""

        #         closest_cols = [[[] for _ in range(6)] for _ in range(closest_n)]

        # transcript_col = []
        # gene_id_col = []
        # symbol_col = []
        # strand_col = []
        # tss_col = []
        # status_col = []

        annotation_cols = [[] for _ in range(6)]

        for _, row in self._df_query.iterrows():
            key = row.name  # (row["Chromosome"], row["Start"], row["End"])
            # convert frozensets back to dict
            annotations = [
                json.loads(x) for x in sorted(annotation_map.get(key, set()))
            ]

            add_annotation_for_location_to_cols(annotations, annotation_cols)

            # if closest_n > 0:
            #     for i in range(1, closest_n + 1):
            #         annotations = [
            #             json.loads(x)
            #             for x in sorted(closest_annotation_map[i].get(key, set()))
            #         ]
            #         add_annotation_for_location_to_cols(
            #             annotations, closest_cols[i - 1]
            #         )

        self._df_query["Transcript Id"] = annotation_cols[0]
        self._df_query["Gene Id"] = annotation_cols[1]
        self._df_query["Gene Symbol"] = annotation_cols[2]
        self._df_query["Strand"] = annotation_cols[3]
        self._df_query["TSS Distance"] = annotation_cols[4]
        self._df_query[self._prom_header] = annotation_cols[5]

        # if closest_n > 0:
        #     for i in range(1, closest_n + 1):
        #         df_query[f"#{i} Transcript Id"] = closest_cols[i - 1][0]
        #         df_query[f"#{i} Gene Id"] = closest_cols[i - 1][1]
        #         df_query[f"#{i} Gene Symbol"] = closest_cols[i - 1][2]
        #         df_query[f"#{i} Strand"] = closest_cols[i - 1][3]
        #         df_query[f"#{i} TSS Distance"] = closest_cols[i - 1][4]
        #         df_query[f"#{i} {prom_header}"] = closest_cols[i - 1][5]

        # query_ranges = pr.PyRanges(df_query)
        # nearest = midpoint_ranges.k_nearest(ALL, k=5, suffix="_nearest", nb_cpu=2)

    def annotate_closest_genes(self, closest_n: int = -1):
        # use default if not specified
        if closest_n == -1:
            closest_n = self._closest_n

        print(f"Finding the {closest_n} closest annotations...")
        # keep track of how many closest are assigned at a location
        used_symbols = collections.defaultdict(dict)
        closest_annotation_map = collections.defaultdict(
            lambda: collections.defaultdict(set)
        )

        print("Processing closest gene annotations...")

        self._cursor.execute(
            NEAREST_GENE_JOIN_QUERY, {"max_distance": self._max_distance}
        )

        rows = []
        for row in self._cursor:
            d = row_to_dict(row)

            location = d["location"]
            gene_symbol = d["gene_symbol"]

            if gene_symbol not in used_symbols[location]:
                if len(used_symbols[location]) < closest_n:
                    used_symbols[location][gene_symbol] = (
                        len(used_symbols[location]) + 1
                    )

            # we keep the first n genes we encounter per location
            closest = used_symbols[location].get(gene_symbol, -1)

            if closest != -1:
                rows.append(d)

        # cursor.execute(DROP_INTRONIC_TABLE_SQL)
        # self._cursor.execute(TEMP_INTRONIC_TABLE_SQL)
        # self._cursor.execute(TEMP_INTRONIC_INDEX_SQL)
        self._cursor.execute(DELETE_INTRONIC_TABLE_SQL)

        queries = []
        for idx, row in enumerate(rows):
            queries.append(
                {
                    "location": row["location"],
                    "chr": row["chr"],
                    "midpoint": row["midpoint"],
                    "row_idx": idx,
                    "transcript_id": row["transcript_id"],
                }
            )

        self._cursor.executemany(
            INSERT_INTRONIC_QUERY,
            queries,
        )

        print("Processing closest intronic annotations...")

        # see if we are in a gene
        self._cursor.execute(IS_INTRONIC_QUERY)

        for c in self._cursor:
            if "intronic" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("intronic")

            if rows[c[0]]["tss_distance"] == -145239 and rows[c[0]]["chr"] == "chr1":
                print("r1", c)
                print("r2", rows[c[0]])

        print("Processing closest exonic annotations...")

        # ok, see if we are in an exon
        self._cursor.execute(IS_EXONIC_QUERY)

        for c in self._cursor:
            if "exonic" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("exonic")

        # for row in rows:
        #    row_to_closest_annotation(row, closest_annotation_map, used_symbols)

        print("Processing closest promoter annotations...")

        self._cursor.execute(
            IS_PROMOTER_QUERY,
            {
                "strand": "+",
                "promoter_lim_1": self._promoter_lim[0],
                "promoter_lim_2": self._promoter_lim[1],
            },
        )

        for c in self._cursor:
            if "promoter" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("promoter")

        self._cursor.execute(
            IS_PROMOTER_QUERY,
            {
                "strand": "-",
                "promoter_lim_1": self._promoter_lim[1],
                "promoter_lim_2": self._promoter_lim[0],
            },
        )

        for c in self._cursor:
            if "promoter" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("promoter")

        for row in rows:
            row_to_closest_annotation(row, closest_annotation_map, used_symbols)

        # for i in range(1, closest_n + 1):
        #     self._df_query[f"#{i} Transcript Id"] = ""
        #     self._df_query[f"#{i} Gene Id"] = ""
        #     self._df_query[f"#{i} Gene Symbol"] = ""
        #     self._df_query[f"#{i} Strand"] = ""
        #     self._df_query[f"#{i} TSS Distance"] = ""
        #     self._df_query[f"#{i} {prom_header}"] = ""

        closest_cols = [[[] for _ in range(6)] for _ in range(closest_n)]

        for _, row in self._df_query.iterrows():
            key = row.name  # (row["Chromosome"], row["Start"], row["End"])

            for i in range(1, closest_n + 1):
                annotations = [
                    json.loads(x)
                    for x in sorted(closest_annotation_map[i].get(key, set()))
                ]
                add_annotation_for_location_to_cols(annotations, closest_cols[i - 1])

        for i in range(1, closest_n + 1):
            self._df_query[f"#{i} Transcript Id"] = closest_cols[i - 1][0]
            self._df_query[f"#{i} Gene Id"] = closest_cols[i - 1][1]
            self._df_query[f"#{i} Gene Symbol"] = closest_cols[i - 1][2]
            self._df_query[f"#{i} Strand"] = closest_cols[i - 1][3]
            self._df_query[f"#{i} TSS Distance"] = closest_cols[i - 1][4]
            self._df_query[f"#{i} {self._prom_header}"] = closest_cols[i - 1][5]


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

        self._cursor.execute(TEMP_QUERY_TABLE_SQL)
        self._cursor.execute(TEMP_INDEX_QUERY_TABLE_REGION_SQL)
        self._cursor.execute(TEMP_INDEX_MID_SQL)

        self._cursor.execute(DROP_INTRONIC_TABLE_SQL)
        self._cursor.execute(TEMP_INTRONIC_TABLE_SQL)
        self._cursor.execute(TEMP_INTRONIC_INDEX_SQL)

    def close(self):
        if self._cursor:
            self._cursor.close()
        if self._conn:
            self._conn.close()

    def annotate_genes(
        self,
    ):
        print(f"Annotating {len(self._df_query)} regions using {self._db}")

        annotation_map = collections.defaultdict(set)

        print("Processing introns...")

        self._cursor.execute(INTRONIC_JOIN_QUERY)

        # for row in cursor:
        #    row_to_annotation(row, annotation_map)

        rows = []
        for row in self._cursor:
            d = row_to_dict(row)

            rows.append(d)

        queries = []
        for idx, row in enumerate(rows):
            queries.append(
                {
                    "row_idx": idx,
                    "location": row["location"],
                    "chr": row["chr"],
                    "midpoint": row["midpoint"],
                    "transcript_id": row["transcript_id"],
                }
            )

        self._cursor.executemany(
            INSERT_INTRONIC_QUERY,
            queries,
        )

        print("Processing exons...")

        # find out which intronic regions are exonic

        self._cursor.execute(IS_EXONIC_QUERY)

        for c in self._cursor:
            if "exonic" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("exonic")

        for row in rows:
            # print(row)
            row_to_annotation(row, annotation_map)

        print("Processing promoters...")

        self._cursor.execute(
            PROMOTER_QUERY,
            {
                "strand": "+",
                "promoter_lim_1": self._promoter_lim[0],
                "promoter_lim_2": self._promoter_lim[1],
            },
        )

        for row in self._cursor:
            row_to_annotation(row_to_dict(row), annotation_map)

        self._cursor.execute(
            PROMOTER_QUERY,
            {
                "strand": "-",
                "promoter_lim_1": self._promoter_lim[1],
                "promoter_lim_2": self._promoter_lim[0],
            },
        )

        for row in self._cursor:
            row_to_annotation(row_to_dict(row), annotation_map)

        print("Adding annotations...")

        for row in self._cursor:
            row_to_annotation(row_to_dict(row), annotation_map)

        # if closest_n > 0:
        #     self.annotate_closest_genes(
        #         cursor, closest_n, max_distance, promoter_lim, closest_annotation_map
        #     )

        # add columns we are going to fill

        # self._df_query["Transcript Id"] = ""
        # self._df_query["Gene Id"] = ""
        # self._df_query["Gene Symbol"] = ""
        # self._df_query["Strand"] = ""
        # self._df_query["TSS Distance"] = ""
        # self._df_query[prom_header] = ""

        # if closest_n > 0:
        #     for i in range(1, closest_n + 1):
        #         df_query[f"#{i} Transcript Id"] = ""
        #         df_query[f"#{i} Gene Id"] = ""
        #         df_query[f"#{i} Gene Symbol"] = ""
        #         df_query[f"#{i} Strand"] = ""
        #         df_query[f"#{i} TSS Distance"] = ""
        #         df_query[f"#{i} {prom_header}"] = ""

        #         closest_cols = [[[] for _ in range(6)] for _ in range(closest_n)]

        # transcript_col = []
        # gene_id_col = []
        # symbol_col = []
        # strand_col = []
        # tss_col = []
        # status_col = []

        annotation_cols = [[] for _ in range(6)]

        for _, row in self._df_query.iterrows():
            key = row.name  # (row["Chromosome"], row["Start"], row["End"])
            # convert frozensets back to dict
            annotations = [
                json.loads(x) for x in sorted(annotation_map.get(key, set()))
            ]

            add_annotation_for_location_to_cols(annotations, annotation_cols)

            # if closest_n > 0:
            #     for i in range(1, closest_n + 1):
            #         annotations = [
            #             json.loads(x)
            #             for x in sorted(closest_annotation_map[i].get(key, set()))
            #         ]
            #         add_annotation_for_location_to_cols(
            #             annotations, closest_cols[i - 1]
            #         )

        self._df_query["Transcript Id"] = annotation_cols[0]
        self._df_query["Gene Id"] = annotation_cols[1]
        self._df_query["Gene Symbol"] = annotation_cols[2]
        self._df_query["Strand"] = annotation_cols[3]
        self._df_query["TSS Distance"] = annotation_cols[4]
        self._df_query[self._prom_header] = annotation_cols[5]

        # if closest_n > 0:
        #     for i in range(1, closest_n + 1):
        #         df_query[f"#{i} Transcript Id"] = closest_cols[i - 1][0]
        #         df_query[f"#{i} Gene Id"] = closest_cols[i - 1][1]
        #         df_query[f"#{i} Gene Symbol"] = closest_cols[i - 1][2]
        #         df_query[f"#{i} Strand"] = closest_cols[i - 1][3]
        #         df_query[f"#{i} TSS Distance"] = closest_cols[i - 1][4]
        #         df_query[f"#{i} {prom_header}"] = closest_cols[i - 1][5]

        # query_ranges = pr.PyRanges(df_query)
        # nearest = midpoint_ranges.k_nearest(ALL, k=5, suffix="_nearest", nb_cpu=2)

    def annotate_closest_genes(self, closest_n: int = -1):
        # use default if not specified
        if closest_n == -1:
            closest_n = self._closest_n

        print(f"Finding the {closest_n} closest annotations...")
        # keep track of how many closest are assigned at a location
        used_symbols = collections.defaultdict(dict)
        closest_annotation_map = collections.defaultdict(
            lambda: collections.defaultdict(set)
        )

        print("Processing closest gene annotations...")

        self._cursor.execute(
            NEAREST_GENE_JOIN_QUERY, {"max_distance": self._max_distance}
        )

        rows = []
        for row in self._cursor:
            d = row_to_dict(row)

            location = d["location"]
            gene_symbol = d["gene_symbol"]

            if gene_symbol not in used_symbols[location]:
                if len(used_symbols[location]) < closest_n:
                    used_symbols[location][gene_symbol] = (
                        len(used_symbols[location]) + 1
                    )

            # we keep the first n genes we encounter per location
            closest = used_symbols[location].get(gene_symbol, -1)

            if closest != -1:
                rows.append(d)

        # cursor.execute(DROP_INTRONIC_TABLE_SQL)
        # self._cursor.execute(TEMP_INTRONIC_TABLE_SQL)
        # self._cursor.execute(TEMP_INTRONIC_INDEX_SQL)
        self._cursor.execute(DELETE_INTRONIC_TABLE_SQL)

        queries = []
        for idx, row in enumerate(rows):
            queries.append(
                {
                    "location": row["location"],
                    "chr": row["chr"],
                    "midpoint": row["midpoint"],
                    "row_idx": idx,
                    "transcript_id": row["transcript_id"],
                }
            )

        self._cursor.executemany(
            INSERT_INTRONIC_QUERY,
            queries,
        )

        print("Processing closest intronic annotations...")

        # see if we are in a gene
        self._cursor.execute(IS_INTRONIC_QUERY)

        for c in self._cursor:
            if "intronic" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("intronic")

            if rows[c[0]]["tss_distance"] == -145239 and rows[c[0]]["chr"] == "chr1":
                print("r1", c)
                print("r2", rows[c[0]])

        print("Processing closest exonic annotations...")

        # ok, see if we are in an exon
        self._cursor.execute(IS_EXONIC_QUERY)

        for c in self._cursor:
            if "exonic" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("exonic")

        # for row in rows:
        #    row_to_closest_annotation(row, closest_annotation_map, used_symbols)

        print("Processing closest promoter annotations...")

        self._cursor.execute(
            IS_PROMOTER_QUERY,
            {
                "strand": "+",
                "promoter_lim_1": self._promoter_lim[0],
                "promoter_lim_2": self._promoter_lim[1],
            },
        )

        for c in self._cursor:
            if "promoter" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("promoter")

        self._cursor.execute(
            IS_PROMOTER_QUERY,
            {
                "strand": "-",
                "promoter_lim_1": self._promoter_lim[1],
                "promoter_lim_2": self._promoter_lim[0],
            },
        )

        for c in self._cursor:
            if "promoter" not in rows[c[0]]["type"]:
                rows[c[0]]["type"].append("promoter")

        for row in rows:
            row_to_closest_annotation(row, closest_annotation_map, used_symbols)

        # for i in range(1, closest_n + 1):
        #     self._df_query[f"#{i} Transcript Id"] = ""
        #     self._df_query[f"#{i} Gene Id"] = ""
        #     self._df_query[f"#{i} Gene Symbol"] = ""
        #     self._df_query[f"#{i} Strand"] = ""
        #     self._df_query[f"#{i} TSS Distance"] = ""
        #     self._df_query[f"#{i} {prom_header}"] = ""

        closest_cols = [[[] for _ in range(6)] for _ in range(closest_n)]

        for _, row in self._df_query.iterrows():
            key = row.name  # (row["Chromosome"], row["Start"], row["End"])

            for i in range(1, closest_n + 1):
                annotations = [
                    json.loads(x)
                    for x in sorted(closest_annotation_map[i].get(key, set()))
                ]
                add_annotation_for_location_to_cols(annotations, closest_cols[i - 1])

        for i in range(1, closest_n + 1):
            self._df_query[f"#{i} Transcript Id"] = closest_cols[i - 1][0]
            self._df_query[f"#{i} Gene Id"] = closest_cols[i - 1][1]
            self._df_query[f"#{i} Gene Symbol"] = closest_cols[i - 1][2]
            self._df_query[f"#{i} Strand"] = closest_cols[i - 1][3]
            self._df_query[f"#{i} TSS Distance"] = closest_cols[i - 1][4]
            self._df_query[f"#{i} {self._prom_header}"] = closest_cols[i - 1][5]

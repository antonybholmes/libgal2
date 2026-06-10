import collections
import itertools
import json
import sqlite3
import sys
from tkinter import LEFT
from typing import Union

import numpy as np
import pandas as pd
import pyranges as pr

from ..genomic.location import NA_DOT

SEP = "|"
PROMOTER_LIM = [2000, 1000]
TOP = 40

# Example: Execute a SELECT query
# TRANSCRIPT_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.symbol,
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

# INTRAGENIC_JOIN_QUERY = f"""
# SELECT DISTINCT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.symbol,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'intronic' AS type,
#     q.midpoint - g.tss AS tss_dist
#     (
#         (g.strand = '+' AND (g.tss - :promoter_lim_1) <= q.midpoint AND (g.tss + :promoter_lim_2) >= q.midpoint)
#         OR
#         (g.strand = '-' AND (g.tss - :promoter_lim_2) <= q.midpoint AND (g.tss + :promoter_lim_1) >= q.midpoint)
#     ) as is_promoter
# FROM query_regions q
# JOIN gtf g ON
#     g.feature = 'transcript' AND
#     g.seqname = q.chr AND
#     (
#         (g.strand = '+' AND (g.tss - :promoter_lim_1) <= q.midpoint AND (g.tss + :promoter_lim_2) >= q.midpoint)
#         OR
#         (g.strand = '-' AND (g.tss - :promoter_lim_2) <= q.midpoint AND (g.tss + :promoter_lim_1) >= q.midpoint)
#     )
# ORDER BY q.location
# """

# INTRAGENIC_JOIN_QUERY = f"""
# INSERT INTO query_closest_genes (location,
#     chr,
#     start,
#     end,
#     midpoint,
#     strand,
#     tss_dist,
#     gene_id,
#     symbol,
#     transcript_id,
#     is_intragenic,
#     is_promoter,
#     gene_rank)
# SELECT q.location,
#     q.chr,
#     q.start,
#     q.end,
#     q.midpoint,
#     g.strand,
#     CASE
#         WHEN g.strand = '-' THEN g.tss - q.midpoint
#         ELSE q.midpoint - g.tss
#     END AS tss_dist,
#     g.gene_id,
#     g.symbol,
#     g.transcript_id,
#     q.midpoint >= g.start AND q.midpoint <= g.end AS is_intragenic,
#     (
#         (g.strand = '+' AND q.midpoint >= (g.tss - :promoter_lim_1) AND q.midpoint <= (g.tss + :promoter_lim_2))
#         OR
#         (g.strand = '-' AND q.midpoint >= (g.tss - :promoter_lim_2) AND q.midpoint <= (g.tss + :promoter_lim_1))
#     ) as is_promoter,
#     1 AS gene_rank
# FROM query_regions q
# JOIN gtf g ON g.feature = 'transcript' AND
# g.seqname = q.chr AND
# (
#     (g.strand = '+' AND q.midpoint >= (g.start - :promoter_lim_1) AND q.midpoint <= g.end)
#     OR
#     (g.strand = '-' AND q.midpoint >= g.start AND q.midpoint <= (g.end + :promoter_lim_1))
# )
# """

# EXON_QUERY = f"""SELECT DISTINCT g.gene_id,
#     g.symbol,
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

# EXON_JOIN_QUERY = f"""
# SELECT DISTINCT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.symbol,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'exonic' AS type,
#     q.midpoint - g.tss AS tss_dist
# FROM query_regions q
# JOIN gtf g ON
#     g.feature = 'exon' AND
#     g.seqname = q.chr AND
#     g.start <= q.midpoint AND
#     g.end >= q.midpoint
# ORDER BY q.location
# """

# IS_INTRAGENIC_JOIN_QUERY = f"""
# SELECT DISTINCT q.row_idx
# FROM intragenic_query_regions q
# JOIN gtf g ON
#     g.feature = 'transcript' AND
#     g.seqname = q.chr AND
#     g.start <= q.midpoint AND
#     g.end >= q.midpoint AND
#     g.transcript_id = q.transcript_id
# ORDER BY q.row_idx
# """

# IS_EXONIC_JOIN_QUERY = f"""
# SELECT DISTINCT q.row_idx
# FROM intragenic_query_regions q
# JOIN gtf g ON
#     g.feature = 'exon' AND
#     g.seqname = q.chr AND
#     g.start <= q.midpoint AND
#     g.end >= q.midpoint AND
#     g.transcript_id = q.transcript_id
# ORDER BY q.row_idx
# """

# IS_PROMOTER_JOIN_QUERY = f"""
# SELECT DISTINCT q.row_idx
# FROM intragenic_query_regions q
# JOIN gtf g ON
#     g.feature = 'transcript' AND
#     g.seqname = q.chr AND
#     (g.tss - :promoter_lim_1) <= q.midpoint AND
#     (g.tss + :promoter_lim_2) >= q.midpoint AND
#     g.strand = :strand AND
#     g.transcript_id = q.transcript_id
# ORDER BY q.row_idx
# """


# PROMOTER_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.symbol,
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

# PROMOTER_QUERY = f"""
# SELECT DISTINCT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.symbol,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'promoter' AS type,
#     q.midpoint - g.tss AS tss_dist
# FROM query_regions q
# JOIN gtf g ON
#     g.feature = 'transcript' AND
#     g.seqname = q.chr AND
#     (
#         (g.strand = '+' AND (g.tss - :promoter_lim_1) <= q.midpoint AND (g.tss + :promoter_lim_2) >= q.midpoint)
#         OR
#         (g.strand = '-' AND (g.tss - :promoter_lim_2) <= q.midpoint AND (g.tss + :promoter_lim_1) >= q.midpoint)
#     )
# ORDER BY q.location
# """


# NEAREST_intragenic_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.symbol,
#     transcript.transcript_id,
#     transcript.tss,
#     gene.strand,
#     ABS(transcript.tss - :midpoint) AS ab_dist,
#     'intronic'
#     FROM gene, transcript
#     WHERE
#     gene.chr = :chromosome AND
#     gene.id = transcript.gene_id
#     ORDER BY ab_dist, gene.symbol
#     LIMIT 50
#     """


# CLOSEST_GENE_JOIN_QUERY = f"""
# SELECT DISTINCT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.symbol,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'intergenic' AS type,
#     g.start,
#     g.end,
#     q.midpoint - g.tss AS tss_dist
# FROM query_regions q
# JOIN gtf g ON
#     g.feature = 'transcript' AND
#     g.seqname = q.chr
# WHERE ABS(tss_dist) < :max_distance
# ORDER BY q.location, ABS(tss_dist), g.symbol
# """

# CLOSEST_GENE_QUERY = f"""
# SELECT DISTINCT
#     g.gene_id,
#     g.symbol,
#     :midpoint - g.tss AS tss_dist
# FROM gtf g
# WHERE g.feature = 'gene' AND g.seqname = :chromosome
# ORDER BY ABS(tss_dist), g.symbol
# LIMIT :limit
# """

# CLOSEST_GENE_GROUP_BY_QUERY = f"""
# INSERT INTO query_closest_genes (location, chr, start, end, midpoint, strand, tss_dist, gene_id, symbol, gene_rank)
# SELECT DISTINCT location, chr, start, end, midpoint, strand, tss_dist, gene_id, symbol, gene_rank
# FROM (
#     SELECT
#         q.location,
#         q.chr,
#         q.start,
#         q.end,
#         q.midpoint,
#         g.strand,
#         g.tss - q.midpoint AS tss_dist,
#         ABS(q.midpoint - g.tss) AS abs_tss_dist,
#         g.gene_id,
#         g.symbol,
#         ROW_NUMBER() OVER (PARTITION BY q.location ORDER BY abs_tss_dist) AS gene_rank
#     FROM query_regions q
#     JOIN gtf g ON g.feature = 'gene' AND g.seqname = q.chr
# )
# WHERE gene_rank <= :limit
# """


INTRAGENIC_GENES_QUERY = f"""
SELECT
    g.id, 
    c.name AS chr,
    g.start, 
    g.end, 
    g.strand, 
    g.gene_id, 
    g.symbol,
    gt.name AS gene_biotype,
    t.transcript_id,
    t.start,
    t.end,
    t.is_canonical,
    t.is_longest,
    ft.name AS feature_type,
    f.start,
    f.end,
    e.exon_id,
    e.exon_number,
    CASE
        WHEN g.strand = '+' THEN :mid - t.start
        ELSE :mid - t.end
    END AS tss_dist,
    ((g.strand = '+') AND (:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p)) OR
		((g.strand = '-') AND (:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p)) 
    AS in_promoter,
    :start <= f.end AND :end >= f.start AS in_exon,
    :start <= t.end AND :end >= t.start AS is_intragenic
FROM genes as g
JOIN chromosomes AS c ON g.chr_id = c.id
JOIN transcripts AS t ON g.id = t.gene_id
JOIN features AS f ON f.transcript_id = t.id
JOIN feature_types AS ft ON f.feature_type_id = ft.id
JOIN exons AS e ON f.exon_id = e.id
JOIN biotypes AS gt ON g.biotype_id = gt.id
WHERE
    c.name = :chr AND 
    (
        (g.strand = '+' AND (:start <= t.end) AND (:end >= t.start - :prom5p)) OR
        (g.strand = '-' AND (:start <= t.end + :prom5p) AND (:end >= t.start))
    ) AND
    LOWER(g.symbol) NOT LIKE 'ensg%'
ORDER BY
    g.gene_id, t.transcript_id, e.exon_number, ABS(tss_dist) ASC
"""

# INSERT INTO query_closest_genes (location,
#     chr,
#     start,
#     end,
#     midpoint,
#     strand,
#     tss_dist,
#     gene_id,
#     symbol,
#     transcript_id,
#     is_intragenic,
#     in_exon,
#     is_promoter,
#     gene_rank)

CLOSEST_GENE_GROUP_BY_QUERY = f"""
-- First, rank transcripts by distance to the query point within each gene and get the n closest transcripts
WITH ranked_transcripts AS (
    SELECT DISTINCT
        g.id as gene_id,
        t.id AS transcript_id,
        ABS(
            CASE 
                WHEN g.strand ='+' THEN :mid - t.start
                ELSE :mid - t.end
            END
        ) as tss_dist,
        -- Rank transcripts within the SAME gene by closest distance
        ROW_NUMBER() OVER (
            PARTITION BY t.gene_id 
            ORDER BY 
                ABS(
                    CASE 
                        WHEN g.strand ='+' THEN :mid - t.start
                        ELSE :mid - t.end
                    END
                ) ASC
        ) AS gene_transcript_rank
    FROM transcripts AS t
    JOIN genes AS g ON t.gene_id = g.id
    JOIN chromosomes AS c ON g.chr_id = c.id
    JOIN biotypes AS gt ON g.biotype_id = gt.id
    WHERE 
        c.name = :chr AND
        -- avoid annotating to genes with an ENSG symbol as these are likely to be less well characterized 
        LOWER(g.symbol) NOT LIKE 'ensg%'
),
closest_transcripts AS (
    SELECT DISTINCT
        ROW_NUMBER() OVER (ORDER BY ABS(tss_dist) ASC) AS rank,
        transcript_id
    FROM 
        ranked_transcripts
    WHERE
        gene_transcript_rank = 1
    ORDER BY 
        ABS(tss_dist) ASC
    LIMIT 
        :limit	
)
-- Then, join the ranked transcripts back to the query points to get the closest gene information for each query point
SELECT DISTINCT
    ct.rank,
    g.id, 
    c.name AS chr,
    g.start, 
    g.end, 
    g.strand, 
    g.gene_id, 
    g.symbol,
    gt.name AS gene_biotype,
    t.transcript_id,
    t.start,
    t.end,
    t.is_canonical,
    t.is_longest,
    CASE
        WHEN g.strand = '+' THEN :mid - t.start
        ELSE :mid - t.end
    END AS tss_dist,
    ((g.strand = '+') AND (:start <= t.start + :prom3p) AND (:end >= t.start - :prom5p)) OR
		((g.strand = '-') AND (:start <= t.end + :prom5p) AND (:end >= t.end - :prom3p)) 
    AS in_promoter,
    :start <= f.end AND :end >= f.start AS in_exon,
    :start <= t.end AND :end >= t.start AS is_intragenic
FROM genes as g
JOIN chromosomes AS c ON g.chr_id = c.id
JOIN transcripts AS t ON g.id = t.gene_id
JOIN closest_transcripts ct ON ct.transcript_id = t.id
JOIN features AS f ON f.transcript_id = t.id
JOIN feature_types AS ft ON f.feature_type_id = ft.id
JOIN exons AS e ON f.exon_id = e.id
JOIN biotypes AS gt ON g.biotype_id = gt.id
ORDER BY 
    ct.rank ASC
"""

SELECT_CLOSEST_GENES_QUERY = f"""SELECT * FROM query_closest_genes"""


QUERY_COUNT_QUERY = f"""
SELECT COUNT(*) AS count
FROM query_regions
"""

CLOSEST_GENE_GROUP_BY_COUNT_QUERY = f"""
SELECT COUNT(*) AS count
FROM query_closest_genes
"""

CLOSEST_TRANSCRIPT_GROUP_BY_COUNT_QUERY = f"""
SELECT COUNT(*) AS count
FROM query_closest_transcripts
"""


# find all occurences of being exonic, but keep one entry per transcript
# for reference
CLOSEST_IS_INTRAGENIC_QUERY = f"""
SELECT DISTINCT
    q.location, 
    q.gene_id
FROM query_closest_genes q
JOIN gtf g ON g.feature = 'transcript' AND 
    g.gene_id = q.gene_id AND 
    g.start <= q.end AND 
    g.end >= q.start
"""

CLOSEST_IS_PROMOTER_QUERY = f"""
SELECT DISTINCT
    q.location, 
    q.gene_id
FROM query_closest_genes q
JOIN gtf g ON g.feature = 'transcript' AND 
    g.gene_id = q.gene_id AND 
    (
        (g.strand = '+' AND (g.tss - :promoter_lim_1) <= q.midpoint AND (g.tss + :promoter_lim_2) >= q.midpoint) 
        OR
        (g.strand = '-' AND (g.tss - :promoter_lim_2) <= q.midpoint AND (g.tss + :promoter_lim_1) >= q.midpoint)
    )
"""

CLOSEST_IS_EXON_QUERY = f"""
SELECT DISTINCT
    q.location, 
    q.gene_id,
    g.transcript_id
FROM query_closest_genes q
JOIN gtf g ON g.feature = 'exon' AND 
    g.gene_id = q.gene_id AND 
    g.start <= q.end AND 
    g.end >= q.start
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


DELETE_QUERY_TABLE_SQL = f"""
    DELETE FROM query_regions
"""

TEMP_INDEX_QUERY_TABLE_REGION_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, start, end)"
)

TEMP_INDEX_QUERY_TABLE_MID_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, midpoint)"
)

TEMP_INDEX_QUERY_TABLE_LOCATION_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_query_location ON query_regions (location)"
)


INSERT_TEMP_QUERY = f"""
    INSERT INTO query_regions (location, chr, start, end, midpoint, strand)
    VALUES (:location, :chr, :start, :end, :midpoint, :strand)
"""


TEMP_CLOSEST_GENE_TABLE_SQL = f"""
    CREATE TEMP TABLE IF NOT EXISTS query_closest_genes (
    location TEXT NOT NULL,
    chr TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    midpoint INTEGER NOT NULL,
    strand TEXT NOT NULL,
    tss_dist INTEGER NOT NULL,
    gene_id TEXT NOT NULL,
    symbol TEXT NOT NULL,
    transcript_id TEXT NOT NULL,
    is_intragenic BOOLEAN NOT NULL,
    in_exon BOOLEAN NOT NULL,
    is_promoter BOOLEAN NOT NULL,
    gene_rank INTEGER NOT NULL
);
"""

# DELETE_TEMP_CLOSEST_GENE_TABLE_SQL = f"""
#     DELETE FROM query_closest_genes
# """


INSERT_TEMP_CLOSEST_GENES_QUERY = f"""
    INSERT INTO query_closest_genes (location, 
    chr, 
    start, 
    end,
    midpoint,
    strand, 
    tss_dist,
    gene_id, 
    symbol, 
    transcript_id,
    gene_rank)
    VALUES (:location, 
    :chr, 
    :start, 
    :end,
    :midpoint,
    :strand, 
    :tss_dist,
    :gene_id, 
    :symbol, 
    :transcript_id,
    :gene_rank)
"""

TEMP_CLOSEST_GENES_INDEX_SQL = (
    f"""CREATE INDEX IF NOT EXISTS idx_closest_genes ON query_closest_genes (gene_id)"""
)

TEMP_CLOSEST_GENES_DELETE_SQL = f"""
    DELETE FROM query_closest_genes
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

    # group by gene id
    gene_id_map = collections.defaultdict(dict)
    for ann in annotations:
        gene_id = ann["gene_id"]

        if gene_id not in gene_id_map:
            gene_id_map[gene_id] = {
                "gene_id": gene_id,
                "symbol": "",
                "labels": set(),
                "tss_dist": {"d": sys.maxsize, "transcript_id": ""},
                "strand": "+",
            }

        # print(ann)

        # entrez_map[ann["gene_id"]]["gene_id"].add(ann["gene_id"])

        gene_id_map[ann["gene_id"]]["symbol"] = ann["symbol"]
        # entrez_map[ann["gene_id"]]["transcript_id"] = ann["transcript_id"]
        gene_id_map[ann["gene_id"]]["labels"].update(ann["labels"])  # (ann["type"])
        # entrez_map[ann["gene_id"]]["tss_dist"] = ann["tss_dist"]

        # we use all annotations for labelling a region, but to reduce
        # noise, we keep the closest refeq by TSS dist per gene and report
        # only that so that we don't have long lists of refseqs for each
        # gene, which are likely redundant
        if abs(ann["tss_dist"]) < abs(gene_id_map[ann["gene_id"]]["tss_dist"]["d"]):
            gene_id_map[ann["gene_id"]]["tss_dist"]["d"] = ann["tss_dist"]
            gene_id_map[ann["gene_id"]]["tss_dist"]["transcript_id"] = ann[
                "transcript_id"
            ]

        gene_id_map[ann["gene_id"]]["strand"] = ann["strand"]

    if len(gene_id_map) > 0:
        gene_id_sorted_annotations = []
        for gene_id in sorted(gene_id_map):
            # collapse tss distance
            ann = gene_id_map[gene_id]

            ann["labels"] = fix_labels(ann["labels"])

            gene_id_sorted_annotations.append(ann)

        # transcript id
        annotation_cols["transcript_id"] = SEP.join(
            ann["tss_dist"]["transcript_id"] for ann in gene_id_sorted_annotations
        )

        # gene id
        annotation_cols["gene_id"] = format_col_str(
            "gene_id", gene_id_sorted_annotations
        )

        annotation_cols["symbol"] = format_col_str("symbol", gene_id_sorted_annotations)

        annotation_cols["strand"] = format_col_str("strand", gene_id_sorted_annotations)
        annotation_cols["tss_distance"] = SEP.join(
            str(ann["tss_dist"]["d"]) for ann in gene_id_sorted_annotations
        )

        annotation_cols["promoter"] = format_col_set(
            "labels", gene_id_sorted_annotations
        )
    else:
        annotation_cols["transcript_id"] = NA_DOT
        annotation_cols["gene_id"] = NA_DOT
        annotation_cols["symbol"] = NA_DOT
        annotation_cols["strand"] = NA_DOT
        annotation_cols["tss_distance"] = NA_DOT
        annotation_cols["promoter"] = "intergenic"


def row_to_dict(row):
    location = row[0]
    chr = row[1]
    midpoint = row[2]
    gene_id = row[3]
    symbol = row[4]
    transcript_id = row[5]
    transcript_tss = row[6]
    gene_strand = row[7]
    transcript_type = row[8]

    if gene_strand == "+":
        dist = midpoint - transcript_tss
    else:
        dist = transcript_tss - midpoint

    annotation = {
        "location": location,
        "chr": chr,
        "midpoint": midpoint,
        "transcript_id": transcript_id,
        "gene_id": gene_id,
        "symbol": symbol,
        "type": [transcript_type],
        "tss_dist": dist,
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
#     symbol = row[4]
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
#         "symbol": symbol,
#         "type": annot_type,
#         "tss_dist": dist,
#         "strand": gene_strand,
#     }

#     return annotation


def row_to_closest_annotation(
    annotation: dict,
    closest_annotation_map: dict[dict[str]],
    used_symbols: dict[dict[int]],
):

    location = annotation["location"]
    symbol = annotation["symbol"]

    # if symbol not in used_symbols[location]:
    #     if len(used_symbols[location]) < closest_n:
    #         used_symbols[location][symbol] = len(used_symbols[location]) + 1

    # we keep the first 5 genes we encounter
    closest = used_symbols[location].get(symbol, -1)

    # if closest == -1:
    # stop
    #    print("stopping")
    #    return True

    # print(symbol, closest, closest_n)

    if closest != -1:
        closest_annotation_map[closest][location].add(
            json.dumps(annotation, sort_keys=True, separators=(",", ":"))
        )


class DataframeAnnotation:
    def __init__(
        self,
        db: str,
        closest_n: int = 5,
        max_distance: int = 2000000,
        promoter_lim: list[int] = [2000, 1000],
        chunk_size: int = 100000,
    ):
        self._db = db
        self._closest_n = closest_n
        self._max_distance = max_distance
        self._promoter_lim = promoter_lim
        self._conn = None
        self._cursor = None
        # self._queries = []
        self._prom_header = f"Relative To Gene (prom=-{promoter_lim[0]/1000}/+{promoter_lim[1]/1000} kb)"
        self._chunk_size = chunk_size

    def _open(self):
        #  close any existing connections
        self.close()

        print(f"Opening database connection to {self._db}")

        self._conn = sqlite3.connect(self._db)

        # Set the row factory to return named tuples
        self._conn.row_factory = sqlite3.Row

        # Create a cursor object
        self._cursor = self._conn.cursor()

    def close(self):
        if self._cursor:
            self._cursor.close()
        if self._conn:
            self._conn.close()

    def annotate_genes(self, query_file: str, out: str):
        self._open()

        print(f"Annotating regions using {self._db}")

        df_iter = pd.read_csv(
            query_file, sep="\t", header=0, chunksize=self._chunk_size
        )
        c = 1

        first = True
        for df_query in df_iter:
            df_query["Transcript Id"] = NA_DOT
            df_query["Gene Id"] = NA_DOT
            df_query["Gene Symbol"] = NA_DOT
            df_query["Strand"] = NA_DOT
            df_query["TSS Distance"] = NA_DOT
            df_query[self._prom_header] = NA_DOT

            for i, query_row in df_query.iterrows():
                if "Start" in query_row:
                    start = query_row["Start"]
                elif "Start_Position" in query_row:
                    start = query_row["Start_Position"]
                else:
                    raise ValueError("Start column not found")

                if "End" in query_row:
                    end = query_row["End"]
                elif "End_Position" in query_row:
                    end = query_row["End_Position"]
                else:
                    raise ValueError("End column not found")

                location = f"{query_row['Chromosome']}:{start}-{end}"

                mid = int((start + end) / 2)

                params = {
                    "location": location,
                    "chr": query_row["Chromosome"],
                    "start": start,
                    "end": end,
                    "mid": mid,
                    "prom5p": self._promoter_lim[0],
                    "prom3p": self._promoter_lim[1],
                }

                # print(params)

                self._cursor.execute(
                    INTRAGENIC_GENES_QUERY,
                    params,
                )

                # get all the annotations for this location
                annotations = []
                for row in self._cursor:
                    r = dict(row)
                    labels = create_labels(row)
                    r["labels"] = list(sorted(labels))
                    annotations.append(r)

                # convert the annotations to columns by appending
                annotation_cols = collections.defaultdict(str)
                add_annotation_for_location_to_cols(annotations, annotation_cols)

                df_query.at[i, "Transcript Id"] = annotation_cols["transcript_id"]
                df_query.at[i, "Gene Id"] = annotation_cols["gene_id"]
                df_query.at[i, "Gene Symbol"] = annotation_cols["symbol"]
                df_query.at[i, "Strand"] = annotation_cols["strand"]
                df_query.at[i, "TSS Distance"] = annotation_cols["tss_distance"]
                df_query.at[i, self._prom_header] = annotation_cols["promoter"]

                if c % 10000 == 0:
                    print(f"Processed {c} queries")
                c += 1

            df_query.to_csv(
                out, sep="\t", header=first, index=False, mode="w" if first else "a"
            )
            first = False

    def annotate_closest_genes(self, query_file: str, out: str, closest_n: int = 1):
        # use default if not specified
        if closest_n == -1:
            closest_n = self._closest_n

        self._open()

        print(f"Processing {closest_n} closest gene annotations...")

        df_iter = pd.read_csv(
            query_file, sep="\t", header=0, chunksize=self._chunk_size
        )

        c = 1
        first = True
        for df_query in df_iter:
            # add the headers for the closest annotations
            for i in range(1, closest_n + 1):
                df_query[f"#{i} Transcript Id"] = NA_DOT
                df_query[f"#{i} Gene Id"] = NA_DOT
                df_query[f"#{i} Gene Symbol"] = NA_DOT
                df_query[f"#{i} Strand"] = NA_DOT
                df_query[f"#{i} TSS Distance"] = NA_DOT
                df_query[f"#{i} {self._prom_header}"] = NA_DOT

            for query_index, query_row in df_query.iterrows():
                if "Start" in query_row:
                    start = query_row["Start"]
                elif "Start_Position" in query_row:
                    start = query_row["Start_Position"]
                else:
                    raise ValueError("Start column not found")

                if "End" in query_row:
                    end = query_row["End"]
                elif "End_Position" in query_row:
                    end = query_row["End_Position"]
                else:
                    raise ValueError("End column not found")

                location = f"{query_row['Chromosome']}:{start}-{end}"

                mid = int((start + end) / 2)

                # now find n closest genes for each row and add to the dataframe
                params = {
                    "location": location,
                    "chr": query_row["Chromosome"],
                    "start": start,
                    "end": end,
                    "mid": mid,
                    "prom5p": self._promoter_lim[0],
                    "prom3p": self._promoter_lim[1],
                    "limit": closest_n,
                }

                # print(params)

                self._cursor.execute(
                    CLOSEST_GENE_GROUP_BY_QUERY,
                    params,
                )

                # for each annotation there might
                # be multiple labels (e.g. exonic, intronic, promoter) so we need to keep track of those and add them to the dataframe as well
                annotation_map = collections.defaultdict(
                    lambda: collections.defaultdict(list)
                )

                for row in self._cursor:
                    rank = row["rank"]
                    annotation_map[rank]["gene_id"].append(row["gene_id"])
                    annotation_map[rank]["transcript_id"].append(row["transcript_id"])
                    annotation_map[rank]["symbol"].append(row["symbol"])
                    annotation_map[rank]["strand"].append(row["strand"])
                    annotation_map[rank]["tss_dist"].append(row["tss_dist"])
                    annotation_map[rank]["labels"].extend(create_labels(row))

                for rank in sorted(annotation_map):
                    rd = annotation_map[rank]

                    df_query.at[query_index, f"#{rank} Transcript Id"] = rd[
                        "transcript_id"
                    ][0]

                    df_query.at[query_index, f"#{rank} Gene Id"] = rd["gene_id"][0]

                    df_query.at[query_index, f"#{rank} Gene Symbol"] = rd["symbol"][0]

                    df_query.at[query_index, f"#{rank} Strand"] = rd["strand"][0]

                    df_query.at[query_index, f"#{rank} TSS Distance"] = rd["tss_dist"][
                        0
                    ]

                    df_query.at[query_index, f"#{rank} {self._prom_header}"] = SEP.join(
                        fix_labels(rd["labels"])
                    )

                if c % 10000 == 0:
                    print(f"Processed {c} queries")
                c += 1

            df_query.to_csv(
                out, sep="\t", header=first, index=False, mode="w" if first else "a"
            )
            first = False


def create_labels(row) -> set[str]:
    labels = set()  # exon_map[id])

    if row["in_exon"]:
        labels.add("exonic")

    if row["in_promoter"]:
        labels.add("promoter")

    if row["is_intragenic"]:
        labels.add("intragenic")

        if not row["in_exon"]:
            labels.add("intronic")

    if len(labels) == 0:
        labels.add("intergenic")

    return labels


def fix_labels(labels: Union[list[str], set[str]]) -> list[str]:
    """Fix the labels for a region based on the following rules:
    - If a region is exonic, it cannot be intronic, so remove the intronic label if it is present.
    - If a region is intronic or exonic, it must be intragenic, so add the intragenic label and remove the intergenic label if it is present.
    - If a region is not intronic or exonic, it must be intergenic, so add the intergenic label and remove the intragenic label if it is present.
    """
    if isinstance(labels, list):
        labels = set(labels)

    if "exonic" in labels:
        labels.discard("intronic")

    if "intronic" in labels or "exonic" in labels:
        labels.add("intragenic")
        labels.discard("intergenic")
    else:
        # ensure that it is marked as not being in a gene
        labels.add("intergenic")

    return list(sorted(labels))

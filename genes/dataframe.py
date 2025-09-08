import collections
import itertools
import json
import sys
from tkinter import LEFT
from typing import Union
import numpy as np
import pandas as pd
import pyranges as pr
from ..genomic.location import NA
import sqlite3


SEP = "|"
PROMOTER_LIM = [2000, 1000]
TOP = 40

# Example: Execute a SELECT query
# TRANSCRIPT_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.gene_name,
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

INTRAGENIC_JOIN_QUERY = f"""
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
    q.midpoint - g.tss AS tss_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND 
    q.midpoint >= g.start AND 
    q.midpoint <= g.end
ORDER BY q.location
"""

# EXON_QUERY = f"""SELECT DISTINCT g.gene_id,
#     g.gene_name,
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
    q.midpoint - g.tss AS tss_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'exon' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint
ORDER BY q.location
"""

IS_INTRAGENIC_JOIN_QUERY = f"""
SELECT DISTINCT q.row_idx
FROM intragenic_query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint AND
    g.transcript_id = q.transcript_id
ORDER BY q.row_idx
"""

IS_EXONIC_JOIN_QUERY = f"""
SELECT DISTINCT q.row_idx
FROM intragenic_query_regions q
JOIN gtf g ON 
    g.feature = 'exon' AND 
    g.seqname = q.chr AND 
    g.start <= q.midpoint AND 
    g.end >= q.midpoint AND
    g.transcript_id = q.transcript_id
ORDER BY q.row_idx
"""

IS_PROMOTER_JOIN_QUERY = f"""
SELECT DISTINCT q.row_idx
FROM intragenic_query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND
    (g.tss - :promoter_lim_1) <= q.midpoint AND 
    (g.tss + :promoter_lim_2) >= q.midpoint AND 
    g.strand = :strand AND
    g.transcript_id = q.transcript_id
ORDER BY q.row_idx
"""


# PROMOTER_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.gene_name,
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
    q.midpoint - g.tss AS tss_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND 
    g.seqname = q.chr AND
    (
        (g.strand = '+' AND (g.tss - :promoter_lim_1) <= q.midpoint AND (g.tss + :promoter_lim_2) >= q.midpoint) 
        OR
        (g.strand = '-' AND (g.tss - :promoter_lim_2) <= q.midpoint AND (g.tss + :promoter_lim_1) >= q.midpoint)
    )
ORDER BY q.location
"""


# NEAREST_intragenic_QUERY = f"""SELECT DISTINCT gene.gene_id,
#     gene.gene_name,
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


CLOSEST_GENE_JOIN_QUERY = f"""
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
    q.midpoint - g.tss AS tss_dist
FROM query_regions q
JOIN gtf g ON 
    g.feature = 'transcript' AND
    g.seqname = q.chr
WHERE ABS(tss_dist) < :max_distance
ORDER BY q.location, ABS(tss_dist), g.gene_name
"""

CLOSEST_GENE_QUERY = f"""
SELECT DISTINCT
    g.gene_id,
    g.gene_name,
    :midpoint - g.tss AS tss_dist
FROM gtf g
WHERE g.feature = 'gene' AND g.seqname = :chromosome
ORDER BY ABS(tss_dist), g.gene_name
LIMIT :limit
"""

CLOSEST_GENE_GROUP_BY_QUERY = f"""
INSERT INTO query_closest_genes (location, chr, start, end, midpoint, gene_id, gene_name, gene_rank)
SELECT DISTINCT location, chr, start, end, midpoint, gene_id, gene_name, gene_rank
FROM (
    SELECT
        q.location,
        q.chr,
        q.start,
        q.end,
        q.midpoint,
        g.gene_id,
        g.gene_name,
        ROW_NUMBER() OVER (PARTITION BY q.location ORDER BY ABS(q.midpoint - g.tss)) AS gene_rank
    FROM query_regions q
    JOIN gtf g ON g.feature = 'gene' AND g.seqname = q.chr
)
WHERE gene_rank <= :limit
"""

# CLOSEST_GENE_GROUP_BY_QUERY = f"""
# WITH ranked_transcripts AS (
#     SELECT DISTINCT
#         q.location,
#         g.gene_id,
#         q.midpoint,
#         ROW_NUMBER() OVER (PARTITION BY q.location ORDER BY ABS(q.midpoint - g.tss)) AS rn
#     FROM query_regions q
#     JOIN gtf g ON g.feature = 'gene' AND g.seqname = q.chr
# )
# SELECT DISTINCT location, gene_id, midpoint
# FROM ranked_transcripts
# WHERE rn <= :limit;
# """

CLOSEST_GENE_GROUP_BY_COUNT_QUERY = f"""
SELECT COUNT(*) AS count
FROM query_closest_genes
"""

CLOSEST_TRANSCRIPT_GROUP_BY_COUNT_QUERY = f"""
SELECT COUNT(*) AS count
FROM query_closest_transcripts
"""

# CLOSEST_TRANSCRIPT_GROUP_BY_QUERY = f"""
# WITH ranked_transcripts AS (
#     SELECT DISTINCT
#         q.location,
#         q.gene_id,
#         g.gene_name,
#         g.transcript_id,
#         g.strand,
#         'closest' AS type,
#         q.midpoint - g.tss AS tss_dist,
#         ABS(q.midpoint - g.tss) AS abs_tss_dist,
#         ROW_NUMBER() OVER (PARTITION BY q.location, q.gene_id ORDER BY ABS(q.midpoint - g.tss)) AS rank
#     FROM query_closest_genes q
#     JOIN gtf g ON g.feature = 'transcript' AND g.gene_id = q.gene_id
#     )
#     SELECT *
#     FROM ranked_transcripts
#     WHERE rank = 1;
# """

INSERT_CLOSEST_TRANSCRIPT_GROUP_BY_QUERY = f"""
INSERT INTO query_closest_transcripts (location, 
    chr, 
    start, 
    end, 
    midpoint, 
    strand, 
    gene_id,
    gene_name,
    transcript_id, 
    tss_dist, 
    abs_tss_dist, 
    is_intragenic,
    is_promoter,
    gene_rank,
    rank)
SELECT DISTINCT location, 
    chr, 
    start, 
    end, 
    midpoint, 
    strand, 
    gene_id,
    gene_name,
    transcript_id, 
    tss_dist, 
    abs_tss_dist, 
    is_intragenic,
    is_promoter,
    gene_rank,
    rank
FROM (
    SELECT q.location,
        q.chr,
        q.start,
        q.end,
        q.midpoint,
        g.strand,
        q.gene_id,
        q.gene_name,
        g.transcript_id,
        q.midpoint - g.tss AS tss_dist,
        ABS(q.midpoint - g.tss) AS abs_tss_dist,
        q.start <= g.end AND q.end >= g.start AS is_intragenic,
        (
            (g.strand = '+' AND (g.tss - :promoter_lim_1) <= q.midpoint AND (g.tss + :promoter_lim_2) >= q.midpoint) 
            OR
            (g.strand = '-' AND (g.tss - :promoter_lim_2) <= q.midpoint AND (g.tss + :promoter_lim_1) >= q.midpoint)
        ) AS is_promoter,
        q.gene_rank,
        ROW_NUMBER() OVER (PARTITION BY q.location, q.gene_id ORDER BY ABS(q.midpoint - g.tss)) AS rank
    FROM query_closest_genes q
    JOIN gtf g ON g.feature = 'transcript' AND 
        g.seqname = q.chr AND 
        g.gene_id = q.gene_id
)
WHERE rank = 1
"""

# find all occurences of being exonic, but keep one entry per transcript
# for reference
SELECT_CLOSEST_TRANSCRIPTS_QUERY = f"""SELECT * FROM query_closest_transcripts;
"""

# find all occurences of being exonic, but keep one entry per transcript
# for reference
SELECT_CLOSEST_EXONS_QUERY = f"""
WITH ranked_exons AS (
    SELECT DISTINCT
        q.chr, 
        q.location, 
        q.gene_id,
        q.gene_name,
        q.transcript_id, 
        q.strand, 
        q.midpoint, 
        q.tss_dist, 
        q.abs_tss_dist, 
        q.is_intragenic,
        q.is_promoter,
        CASE 
            WHEN g.start IS NOT NULL THEN 1 
            ELSE 0 
        END AS is_exonic,
        ROW_NUMBER() OVER (PARTITION BY q.transcript_id ORDER BY g.exon_number) AS rank
    FROM query_closest_transcripts q
    LEFT JOIN gtf g ON g.feature = 'exon' AND 
        g.gene_id = q.gene_id AND 
        g.transcript_id = q.transcript_id AND 
        q.start <= g.end AND 
        q.end >= g.start
)
SELECT *
FROM ranked_exons
WHERE rank = 1;
"""

COUNT_CLOSEST_EXONS_QUERY = f"""
WITH ranked_exons AS (
    SELECT DISTINCT
        q.chr, 
        q.location, 
        q.gene_id,
        q.gene_name,
        q.transcript_id, 
        q.strand, 
        q.midpoint, 
        q.tss_dist, 
        q.abs_tss_dist, 
        q.is_intragenic,
        q.is_promoter,
        CASE 
            WHEN g.start IS NOT NULL THEN 1 
            ELSE 0 
        END AS is_exonic,
        ROW_NUMBER() OVER (PARTITION BY q.transcript_id ORDER BY g.exon_number) AS rank
    FROM query_closest_transcripts q
    LEFT JOIN gtf g ON g.feature = 'exon' AND 
        g.gene_id = q.gene_id AND 
        g.transcript_id = q.transcript_id AND 
        q.start <= g.end AND 
        q.end >= g.start
)
SELECT COUNT(*)
FROM ranked_exons
WHERE rank = 1;
"""

# IS_EXONIC_QUERY = f"""
# SELECT DISTINCT COUNT(id) as count
# FROM gtf g
# WHERE
#     g.feature = 'exon' AND
#     g.transcript_id = :transcript_id AND
#     g.seqname = :chr AND
#     g.start <= :end AND
#     g.end >= :start
# """

# IS_PROMOTER_QUERY = f"""
# SELECT DISTINCT COUNT(id) as count
# FROM gtf g
# WHERE
#     g.feature = 'transcript' AND
#     g.transcript_id = :transcript_id AND
#     g.seqname = :chr AND
#     (
#         (g.strand = '+' AND (g.tss - :promoter_lim_1) <= :midpoint AND (g.tss + :promoter_lim_2) >= :midpoint)
#         OR
#         (g.strand = '-' AND (g.tss - :promoter_lim_2) <= :midpoint AND (g.tss + :promoter_lim_1) >= :midpoint)
#     );
# """

# IS_INTRA_PROM_EXON_JOIN_QUERY = f"""
# SELECT DISTINCT
#     q.row_idx
#     q.location,
#     t.transcript_id,
#     t.start <= :end AND t.end >= :start as is_intragenic,
#     e.start <= :end AND e.end >= :start AS is_exonic,
#     (
#         (t.strand = '+' AND (t.tss - :promoter_lim_1) <= :midpoint AND (t.tss + :promoter_lim_2) >= :midpoint)
#         OR
#         (t.strand = '-' AND (t.tss - :promoter_lim_2) <= :midpoint AND (t.tss + :promoter_lim_1) >= :midpoint)
#     ) AS is_promoter,
# FROM intragenic_query_regions q
# JOIN gtf t ON
#     t.feature = 'transcript' AND
#     t.transcript_id = q.transcript_id
# JOIN gtf e ON
#     e.feature = 'exon' AND
#     e.transcript_id = t.transcript_id
# ORDER BY q.location;
# """

# NEAREST_intragenic_JOIN_QUERY = f"""
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

# NEAREST_EXON_JOIN_QUERY = f"""
# SELECT DISTINCT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.gene_name,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'exonic' AS type,
#     ABS(q.midpoint - g.tss) as ab_dist
# FROM query_regions q
# JOIN gtf g ON
#     g.feature = 'exon' AND
#     g.seqname = q.chr AND
#     g.start <= q.midpoint AND
#     g.end >= q.midpoint
# ORDER BY q.location, ab_dist, g.gene_name;
# """

# NEAREST_PROMOTER_JOIN_QUERY = f"""
# SELECT DISTINCT
#     q.location,
#     q.chr,
#     q.midpoint,
#     g.gene_id,
#     g.gene_name,
#     g.transcript_id,
#     g.tss,
#     g.strand,
#     'promoter' AS type,
#     g.start,
#     g.end,
#     ABS(q.midpoint - g.tss) as ab_dist
# FROM query_regions q
# JOIN gtf g ON
#     g.feature = 'transcript' AND
#     g.seqname = q.chr AND
#     (g.tss - :promoter_lim_1) <= q.midpoint AND
#     (g.tss + :promoter_lim_2) >= q.midpoint AND
#     g.strand = :strand
# ORDER BY q.location, g.gene_name;
# """


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

TEMP_INTRAGENIC_TABLE_SQL = f"""
    CREATE TEMP TABLE IF NOT EXISTS intragenic_query_regions (
    row_idx INTEGER NOT NULL,
    location TEXT NOT NULL,
    chr TEXT NOT NULL,
    midpoint INTEGER NOT NULL,
    transcript_id TEXT NOT NULL
);
"""


INSERT_INTRAGENIC_QUERY = f"""
    INSERT INTO intragenic_query_regions (row_idx, location, chr, midpoint, transcript_id)
    VALUES (:row_idx, :location, :chr, :midpoint, :transcript_id)
"""

DROP_INTRAGENIC_TABLE_SQL = f"""
    DROP TABLE IF EXISTS intragenic_query_regions
"""

DELETE_INTRAGENIC_TABLE_SQL = f"""
    DELETE FROM intragenic_query_regions
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

# TEMP_intragenic_TABLE_SQL = (
#     "CREATE INDEX IF NOT EXISTS idx_query_regions ON query_regions (chr, start, end)"
# )

# TEMP_LOCATION_INDEX_TABLE_SQL = (
#     "CREATE INDEX IF NOT EXISTS idx_query_regions_location ON query_regions (location)"
# )

TEMP_QUERY_INDEX_SQL = f"""CREATE INDEX IF NOT EXISTS idx_intragenic_query_regions ON intragenic_query_regions (chr, midpoint)"""

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
    gene_id TEXT NOT NULL,
    gene_name TEXT NOT NULL,
    gene_rank INTEGER NOT NULL
);
"""

# DELETE_TEMP_CLOSEST_GENE_TABLE_SQL = f"""
#     DELETE FROM query_closest_genes
# """


# INSERT_TEMP_CLOSEST_GENES_QUERY = f"""
#     INSERT INTO query_closest_genes (location, gene_id, midpoint)
#     VALUES (:location, :gene_id, :midpoint)
# """

TEMP_CLOSEST_GENES_INDEX_SQL = (
    f"""CREATE INDEX IF NOT EXISTS idx_closest_genes ON query_closest_genes (gene_id)"""
)


TEMP_CLOSEST_TRANSCRIPT_TABLE_SQL = f"""
    CREATE TEMP TABLE IF NOT EXISTS query_closest_transcripts (
    location TEXT NOT NULL,
    chr TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    midpoint INTEGER NOT NULL,
    strand TEXT NOT NULL,
    gene_id TEXT NOT NULL,
    gene_name TEXT NOT NULL,
    transcript_id TEXT NOT NULL,
    tss_dist INTEGER NOT NULL,
    abs_tss_dist INTEGER NOT NULL,
    is_intragenic INTEGER NOT NULL,
    is_promoter INTEGER NOT NULL,
    gene_rank INTEGER NOT NULL,
    rank INTEGER NOT NULL
);
"""

DELETE_TEMP_CLOSEST_TRANSCRIPT_TABLE_SQL = f"""
    DELETE FROM query_closest_transcripts
"""


TEMP_CLOSEST_TRANSCRIPT_INDEX_SQL = f"""CREATE INDEX IF NOT EXISTS idx_closest_transcripts ON query_closest_transcripts (gene_id, transcript_id)"""


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
                "gene_name": "",
                "labels": set(),
                "tss_dist": {"d": sys.maxsize, "transcript_id": ""},
                "strand": "+",
            }

        print(ann)

        # entrez_map[ann["gene_id"]]["gene_id"].add(ann["gene_id"])

        gene_id_map[ann["gene_id"]]["gene_name"] = ann["gene_name"]
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

            if "exonic" in ann["labels"]:
                ann["labels"].discard("intronic")

            # if ann["labels"] contains intronic or exonic
            if "intronic" in ann["labels"] or "exonic" in ann["labels"]:
                ann["labels"].add("intragenic")
                ann["labels"].discard("intergenic")
            else:
                # ensure that it is marked as not being in a gene
                ann["labels"].add("intergenic")

            gene_id_sorted_annotations.append(gene_id_map[gene_id])

        # transcript id
        annotation_cols[0].append(
            SEP.join(
                ann["tss_dist"]["transcript_id"] for ann in gene_id_sorted_annotations
            )
        )

        # gene id
        annotation_cols[1].append(format_col_str("gene_id", gene_id_sorted_annotations))

        annotation_cols[2].append(
            format_col_str("gene_name", gene_id_sorted_annotations)
        )

        annotation_cols[3].append(format_col_str("strand", gene_id_sorted_annotations))
        annotation_cols[4].append(
            SEP.join(str(ann["tss_dist"]["d"]) for ann in gene_id_sorted_annotations)
        )

        annotation_cols[5].append(format_col_set("labels", gene_id_sorted_annotations))
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
    gene_name = row[4]
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
        "gene_name": gene_name,
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
#     gene_name = row[4]
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
#         "gene_name": gene_name,
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
    gene_name = annotation["gene_name"]

    # if gene_name not in used_symbols[location]:
    #     if len(used_symbols[location]) < closest_n:
    #         used_symbols[location][gene_name] = len(used_symbols[location]) + 1

    # we keep the first 5 genes we encounter
    closest = used_symbols[location].get(gene_name, -1)

    # if closest == -1:
    # stop
    #    print("stopping")
    #    return True

    # print(gene_name, closest, closest_n)

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

        # Set the row factory to return named tuples
        self._conn.row_factory = sqlite3.Row

        # Create a cursor object
        self._cursor = self._conn.cursor()

        self._cursor.execute(TEMP_QUERY_TABLE_SQL)
        self._cursor.execute(TEMP_INDEX_QUERY_TABLE_REGION_SQL)
        self._cursor.execute(TEMP_INDEX_QUERY_TABLE_MID_SQL)
        self._cursor.execute(TEMP_INDEX_QUERY_TABLE_LOCATION_SQL)

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

        self._cursor.execute(DROP_INTRAGENIC_TABLE_SQL)
        self._cursor.execute(TEMP_INTRAGENIC_TABLE_SQL)
        self._cursor.execute(TEMP_QUERY_INDEX_SQL)

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

        self._cursor.execute(INTRAGENIC_JOIN_QUERY)

        # for row in cursor:
        #    row_to_annotation(row, annotation_map)

        rows = [dict(row) for row in self._cursor]

        queries = []
        for idx, row in enumerate(rows):
            row["labels"] = []
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
            INSERT_INTRAGENIC_QUERY,
            queries,
        )

        print("Processing exons...")

        # find out which intronic regions are exonic

        self._cursor.execute(IS_EXONIC_JOIN_QUERY)

        for c in self._cursor:
            if "exonic" not in rows[c[0]]["labels"]:
                rows[c[0]]["labels"].append("exonic")

        for row in rows:
            # print(row)
            row_to_annotation(row, annotation_map)

        print("Processing promoters...")

        self._cursor.execute(
            PROMOTER_QUERY,
            {
                "promoter_lim_1": self._promoter_lim[0],
                "promoter_lim_2": self._promoter_lim[1],
            },
        )

        for row in self._cursor:
            r = dict(row)
            r["labels"] = []
            row_to_annotation(r, annotation_map)

        # self._cursor.execute(
        #     PROMOTER_QUERY,
        #     {
        #         "strand": "-",
        #         "promoter_lim_1": self._promoter_lim[1],
        #         "promoter_lim_2": self._promoter_lim[0],
        #     },
        # )

        # for row in self._cursor:
        #     row_to_annotation(dict(row), annotation_map)

        print("Adding annotations...")

        for row in self._cursor:
            row_to_annotation(row, annotation_map)

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

        self._cursor.execute(TEMP_CLOSEST_GENE_TABLE_SQL)
        self._cursor.execute(TEMP_CLOSEST_GENES_INDEX_SQL)
        # self._cursor.execute(DELETE_TEMP_CLOSEST_GENE_TABLE_SQL)
        self._cursor.execute(TEMP_CLOSEST_TRANSCRIPT_TABLE_SQL)
        self._cursor.execute(TEMP_CLOSEST_TRANSCRIPT_INDEX_SQL)

        print(f"Finding the {closest_n} closest annotations...")
        # keep track of how many closest are assigned at a location
        used_symbols = collections.defaultdict(dict)

        closest_annotation_map = collections.defaultdict(
            lambda: collections.defaultdict(set)
        )

        print(f"Processing {closest_n} closest gene annotations...")

        self._cursor.execute(
            CLOSEST_GENE_GROUP_BY_QUERY,
            {
                "limit": closest_n,
            },
        )

        self._cursor.execute(
            CLOSEST_GENE_GROUP_BY_COUNT_QUERY,
        )

        count = self._cursor.fetchone()["count"]
        print(f"Found {count} closest gene annotations")

        print(f"Processing {closest_n} closest transcript annotations...")

        self._cursor.execute(
            INSERT_CLOSEST_TRANSCRIPT_GROUP_BY_QUERY,
            {
                "promoter_lim_1": self._promoter_lim[0],
                "promoter_lim_2": self._promoter_lim[1],
            },
        )

        self._cursor.execute(
            CLOSEST_TRANSCRIPT_GROUP_BY_COUNT_QUERY,
        )

        count = self._cursor.fetchone()["count"]
        print(f"Found {count} closest transcript annotations")

        # print(f"Processing {closest_n} closest exon annotations...")

        # self._cursor.execute(
        #     COUNT_CLOSEST_EXONS_QUERY,
        # )

        # count = self._cursor.fetchone()["count"]
        # print(f"Found {count} closest exon annotations")

        self._cursor.execute(
            SELECT_CLOSEST_TRANSCRIPTS_QUERY,
        )

        closest_genes = collections.defaultdict(list)

        current_location = None

        # self._df_query[f"#{i} Transcript Id"] = closest_cols[i - 1][0]
        #     self._df_query[f"#{i} Gene Id"] = closest_cols[i - 1][1]
        #     self._df_query[f"#{i} Gene Symbol"] = closest_cols[i - 1][2]
        #     self._df_query[f"#{i} Strand"] = closest_cols[i - 1][3]
        #     self._df_query[f"#{i} TSS Distance"] = closest_cols[i - 1][4]
        #     self._df_query[f"#{i} {self._prom_header}"] = closest_cols[i - 1][5]

        # transcript_ids = []
        # gene_ids = []
        # gene_names = []
        # strands = []
        # tss_dists = []
        # labels = []
        data = np.full(
            (self._df_query.shape[0], 1 + closest_n * 6), "n/a", dtype=object
        )

        # for idx, row in enumerate(self._cursor):
        #    pass

        # print(idx)

        print("Writing closest annotations...")

        for idx, row in enumerate(self._cursor):
            d = dict(row)

            # print(d)

            row_idx = int(idx / closest_n)

            if row_idx >= data.shape[0]:
                print(
                    f"Warning: more closest annotations than query regions {data.shape[0]} < {row_idx} {idx}"
                )
                break

            types = set()

            # if d["is_intragenic"]:
            #     types.add("intragenic")
            # if d["is_exonic"]:
            #     types.add("exonic")
            # if d["is_promoter"]:
            #     types.add("promoter")

            if "intragenic" in types and "exonic" not in types:
                types.add("intronic")

            if len(types) == 0:
                types.add("intergenic")

            closest_idx = d["gene_rank"] - 1  # idx % closest_n
            closest_block = closest_idx * 6 + 1

            data[row_idx, 0] = d["location"]

            data[row_idx, closest_block] = d["transcript_id"]
            data[row_idx, closest_block + 1] = d["gene_id"]
            data[row_idx, closest_block + 2] = d["gene_name"]
            data[row_idx, closest_block + 3] = d["strand"]
            data[row_idx, closest_block + 4] = d["tss_dist"]
            data[row_idx, closest_block + 5] = ",".join(sorted(types))

            if idx % 1000 == 0:
                print(f"Processed {idx} closest annotations...")

        #     if current_location is None:
        #         current_location = d["location"]

        #     if current_location != d["location"]:
        #         data.append(
        #             [d["location"]]
        #             + list(
        #                 itertools.chain(
        #                     [
        #                         [
        #                             transcript_ids[i],
        #                             gene_ids[i],
        #                             gene_names[i],
        #                             strands[i],
        #                             tss_dists[i],
        #                             labels[i],
        #                         ]
        #                         for i in range(closest_n)
        #                     ]
        #                 )
        #             )
        #         )

        #         current_location = d["location"]

        #         transcript_ids = []
        #         gene_ids = []
        #         gene_names = []
        #         strands = []
        #         tss_dists = []
        #         labels = []
        #     else:
        #         transcript_ids.append(d["transcript_id"])
        #         gene_ids.append(d["gene_id"])
        #         gene_names.append(d["gene_name"])
        #         strands.append(d["strand"])
        #         tss_dists.append(d["tss_dist"])
        #         labels.append(d["type"])

        #     closest_genes[d["location"]].append(d)

        # # add the last entry
        # data.append(
        #     [d["location"]]
        #     + list(
        #         itertools.chain(
        #             [
        #                 [
        #                     transcript_ids[i],
        #                     gene_ids[i],
        #                     gene_names[i],
        #                     strands[i],
        #                     tss_dists[i],
        #                     labels[i],
        #                 ]
        #                 for i in range(closest_n)
        #             ]
        #         )
        #     )
        # )

        print(data.shape)

        headers = ["Location"]

        for i in range(1, closest_n + 1):
            headers.append(f"#{i} Transcript Id")
            headers.append(f"#{i} Gene Id")
            headers.append(f"#{i} Gene Symbol")
            headers.append(f"#{i} Strand")
            headers.append(f"#{i} TSS Distance")
            headers.append(f"#{i} {self._prom_header}")

        closest_df = pd.DataFrame(
            data,
            columns=headers,
        )

        closest_df = closest_df.set_index("Location")
        closest_df.index.name = "Location"

        print(closest_df.shape, self._df_query.shape)

        print(closest_df.index)

        self._df_query = pd.concat([self._df_query, closest_df], axis=1)

        # self._cursor.executemany(
        #     INSERT_TEMP_CLOSEST_GENES_QUERY,
        #     closest_genes,
        # )

        # print(len(self._queries), len(closest_genes))

        # self._cursor.execute(
        #     CLOSEST_TRANSCRIPT_GROUP_BY_QUERY,
        # )

        # closest_genes = []

        # for row in self._cursor:
        #     closest_genes.append(dict(row))

        print(len(self._queries), len(closest_genes))

        # queries_closest_genes = []

        # for qi, query in enumerate(self._queries):
        #     closest_genes = []

        #     self._cursor.execute(
        #         CLOSEST_GENE_QUERY,
        #         {
        #             "midpoint": query["midpoint"],
        #             "chromosome": query["chr"],
        #             "limit": closest_n,
        #         },
        #     )

        #     rows = [dict(row) for row in self._cursor]

        #     # print(query, len(rows))

        #     for closest_gene in rows:
        #         self._cursor.execute(
        #             CLOSEST_TRANSCRIPTS_QUERY,
        #             {
        #                 "gene_id": closest_gene["gene_id"],
        #                 "midpoint": query["midpoint"],
        #                 "chr": query["chr"],
        #                 "start": query["start"],
        #                 "end": query["end"],
        #             },
        #         )

        #         row = dict(self._cursor.fetchone())
        #         row["labels"] = set()

        #         if row["is_intragenic"]:
        #             row["labels"].add("intragenic")

        #         # see if exonic
        #         self._cursor.execute(
        #             IS_EXONIC_QUERY,
        #             {
        #                 "transcript_id": row["transcript_id"],
        #                 "midpoint": query["midpoint"],
        #                 "chr": query["chr"],
        #                 "start": query["start"],
        #                 "end": query["end"],
        #             },
        #         )
        #         is_exonic = self._cursor.fetchone()["count"] > 0
        #         if is_exonic:
        #             row["labels"].add("exonic")

        #         self._cursor.execute(
        #             IS_PROMOTER_QUERY,
        #             {
        #                 "transcript_id": row["transcript_id"],
        #                 "midpoint": query["midpoint"],
        #                 "chr": query["chr"],
        #                 "start": query["start"],
        #                 "end": query["end"],
        #                 "promoter_lim_1": self._promoter_lim[0],
        #                 "promoter_lim_2": self._promoter_lim[1],
        #             },
        #         )
        #         is_promoter = self._cursor.fetchone()["count"] > 0
        #         if is_promoter and "promoter" not in row["labels"]:
        #             row["labels"].add("promoter")

        #         if "intragenic" in row["labels"] and "exonic" not in row["labels"]:
        #             row["labels"].add("intronic")

        #         closest_genes.append(row)

        #     print(len(closest_genes))
        #     queries_closest_genes.append(closest_genes)

        #     if (qi + 1) % 1000 == 0:
        #         print(f"Processed {qi + 1} queries...")

        # closest_cols = [[[] for _ in range(6)] for _ in range(closest_n)]

        # for qi, row in self._df_query.iterrows():
        #     key = row.name  # (row["Chromosome"], row["Start"], row["End"])

        #     for i in range(1, closest_n + 1):
        #         annotations = [
        #             json.loads(x)
        #             for x in sorted(closest_annotation_map[i].get(key, set()))
        #         ]
        #         # add_annotation_for_location_to_cols(annotations, closest_cols[i - 1])
        #         annotations = closest_genes[qi]

        #         # transcript id
        #         closest_cols[i - 1][0].append(
        #             SEP.join([ann["transcript_id"] for ann in annotations])
        #         )

        #         closest_cols[i - 1][1].append(
        #             SEP.join([ann["gene_id"] for ann in annotations])
        #         )

        #         closest_cols[i - 1][2].append(
        #             SEP.join([ann["gene_name"] for ann in annotations])
        #         )

        #         closest_cols[i - 1][3].append(
        #             SEP.join([ann["strand"] for ann in annotations])
        #         )

        #         closest_cols[i - 1][4].append(
        #             SEP.join([str(ann["tss_dist"]) for ann in annotations])
        #         )

        #         closest_cols[i - 1][5].append(
        #             SEP.join([",".join(ann["labels"]) for ann in annotations])
        #         )

        # for i in range(1, closest_n + 1):
        #     self._df_query[f"#{i} Transcript Id"] = closest_cols[i - 1][0]
        #     self._df_query[f"#{i} Gene Id"] = closest_cols[i - 1][1]
        #     self._df_query[f"#{i} Gene Symbol"] = closest_cols[i - 1][2]
        #     self._df_query[f"#{i} Strand"] = closest_cols[i - 1][3]
        #     self._df_query[f"#{i} TSS Distance"] = closest_cols[i - 1][4]
        #     self._df_query[f"#{i} {self._prom_header}"] = closest_cols[i - 1][5]

        return self._df_query

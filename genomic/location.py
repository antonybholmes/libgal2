# -*- coding: utf-8 -*-
"""
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

Copyright (C) 2022 Antony Holmes.
"""

from abc import ABC, abstractmethod
import collections
from enum import Enum
from functools import total_ordering
from importlib.metadata import metadata
import re
import sys
from typing import Any, Iterable, Mapping, Union


from .. import species

TELOMERE_SIZE = 100000
NA = "n/a"


def one_base_loc(start: int, end: int) -> tuple[int, int]:
    """
    Ensures start and end are 1 based and end is >= start

    Args:
            start (int): start
            end (int): end

    Returns:
            tuple[int, int]: modified start and end
    """
    start = max(1, start)
    # end cannot be the start
    end = max(start, end)
    return (start, end)


def promoter_type(prom_ext_5p: int = 2000, prom_ext_3p: int = 1000):
    return f"prom=-{prom_ext_5p / 1000}/+{prom_ext_3p / 1000}kb"


def promoter_heading(prom_ext_5p: int = 2000, prom_ext_3p: int = 1000):
    return f"({promoter_type(prom_ext_5p, prom_ext_3p)})"


def location_string(chr: str, start: int, end: int, strand: str = "+") -> str:
    """
    Returns a standardized string representation of a genomic location.

    Args:
            chr (str):    chromosome
            start (int):  1-based start location
            end (int):    1-based end location
    Returns:
            str: A chr:start-end formatted string.
    """
    start, end = one_base_loc(start, end)

    return f"{chr}:{str(start)}-{str(end)}"

    # if strand == "+":
    #     return f"{chr}:{str(start)}-{str(end)}"
    # else:
    #     return f"{chr}:{str(end)}-{str(start)}"


# class Strand(Enum):
#     PLUS = 1
#     MINUS = 2


def chr_to_index(chr: str) -> int:
    """
    Converts chr1,chr2... into a numerical value e.g. chr1->1,chr2->2 etc.
    x,y,m are assigned 1000,2000,3000 respectively to preserve order and allow
    for species with differing numbers of chromosomes (assuming < 1000), but this
    easily covers human, mouse and other common mammals.
    """
    chr = (
        chr.lower()
        .replace("chr", "")
        .replace("x", "1000")
        .replace("y", "2000")
        .replace("m", "3000")
    )

    if chr.isdigit():
        return int(chr)
    else:
        return -1


@total_ordering
class Location:
    """
    Represents a genomic location.
    """

    def __init__(self, chr: str, start: int, end: int, strand: str = "+"):
        """
        Creates a new location object

        Args:
                chr (str): 1 based chromosome, e.g. "chr2"
                start (int): 1 based start location.
                end (int): 1 based end location.
        """
        self._chr = chr  # .lower()

        # turn chr into a numerical id for easier sorting
        self._index = chr_to_index(chr)

        self._start, self._end = one_base_loc(start, end)
        self._length = self._end - self._start + 1
        self._strand = strand

    @property
    def chr(self) -> str:
        return self._chr

    @property
    def index(self) -> int:
        return self._index

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    @property
    def mid(self) -> int:
        return int((self._start + self._end) / 2)

    @property
    def length(self) -> int:
        return self._length

    @property
    def strand(self) -> str:
        return self._strand

    def to_string(self) -> str:
        return self.__str__()

    @property
    def values(self) -> list[str, int, int]:
        return [self.chr, self.start, self.end]

    @property
    def bed(self) -> list[str, int, int]:
        return [self.chr, self.start - 1, self.end]

    @property
    def bed_str(self) -> str:
        return "\t".join([str(x) for x in self.bed])

    @property
    def id(self) -> str:
        return "-".join([str(x) for x in self.values])

    def __str__(self) -> str:
        # , self.strand)
        return location_string(self.chr, self.start, self.end)

    def __repr__(self) -> str:
        return self.__str__()

    def __len__(self) -> int:
        return self._length

    def __hash__(self):
        return hash(self.__str__())

    def __eq__(self, other):

        return (
            self.chr == other.chr
            and self.start == other.start
            and self.end == other.end
        )

    def __lt__(self, other):
        if not isinstance(other, Location):
            return True

        if self.index != other.index:
            return self.index < other.index

        # chr same so test positions

        if self.start != other.start:
            return self.start < other.start

        return self.end < other.end


def sort_locations(locations: Iterable[Location]):
    """Sort locations by chr and start

    Args:
        locations (Iterable[Location]): _description_

    Returns:
        _type_: _description_
    """

    m = collections.defaultdict(lambda: collections.defaultdict(list[Location]))

    for l in locations:
        m[l.chr][l.start].append(l)

    ret = []

    for chr in sorted(m):
        for start in sorted(m[chr]):
            ret.extend(m[chr][start])

    return ret


def mid(location: Location):
    if location is None:
        return -1

    return (location.start + location.end) / 2


def mid_point(location: Location) -> int:
    return int(mid(location))


def center_location(location: Location) -> Location:
    """Returns the mid point of a location as a location object.

    Args:
            location (Location): Genomic location

    Returns:
            Location:  Genomic location
    """

    c = mid_point(location)

    return Location(location.chr, c, c)


def mid_point_dist(l1: Location, l2: Location) -> int:
    if l1 is None or l2 is None:
        return -1

    return abs(mid_point(l1) - mid_point(l2))


def gap_dist(l1: Location, l2: Location) -> int:
    """Closest distance between max start and min end

    Args:
            l1 (Location): _description_
            l2 (Location): _description_

    Returns:
            int: gap distance
    """
    if l1 is None or l2 is None:
        return -1

    return min(
        mid_point_dist(l1, l2), abs(max(l1.start, l2.start) - min(l1.end, l2.end))
    )


class Feature(Location):
    """
    Represents a named genomic location with multiple ids.
    """

    def __init__(self, chr: str, start: int, end: int, strand: str = "+"):
        """_summary_

        Args:
                chr (str): chromosome
                start (int): 1 based start location
                end (int): 1 based end location
        """
        super().__init__(chr, start, end, strand)
        self._id_map = {}

    def add_id(self, id: str, name: str):
        self._id_map[id] = name

    def get_id(self, name: str) -> str:
        """
        Returns a named id associated with the gene or n/a if not present

        Args:
                name:   name of id, e.g. "refseq"
        Returns:
                Named argument.
        """
        return self._id_map.get(name, NA)


class Chromosomes:
    """
    Chromosome sizes
    """

    def __init__(self, file: str):
        self._sizes = collections.defaultdict(int)

        f = open(file, "r")

        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            tokens = line.split("\t")

            chr = tokens[0]
            size = int(tokens[1])

            self._sizes[chr] = size

        f.close()

    def get_size(self, chr):
        return self._sizes[chr]


class FeatureSet:
    def __init__(self, start):
        self._start: int = start
        self._values: list[Location] = []

    @property
    def start(self):
        return self._start

    @property
    def values(self) -> list[Location]:
        return self._values

    def add(self, value: Location):
        self._values.append(value)

    def __iter__(self):
        return iter(self._values)


class GenomicSearch(ABC):
    def get_closest_feature(self, location: Location) -> Location:
        featureSet: FeatureSet = self.get_closest_features(location)

        if featureSet is None:
            print(location, featureSet)
            return None

        min_d = sys.maxsize
        ret: Location = None

        for v in featureSet:
            d = abs(mid(location) - mid(v))

            if d < min_d:
                ret = v
                min_d = d

        return ret

    def get_closest_features(self, location: Location) -> FeatureSet:
        featureSets: list[FeatureSet] = self.get_features(location)

        min_d = sys.maxsize
        ret = None

        for featureSet in featureSets:
            d = abs(mid(location) - featureSet.start)

            if d < min_d:
                ret = featureSet
                min_d = d

        return ret

    @abstractmethod
    def get_features(self, location: Location) -> list[FeatureSet]: ...


class GappedSearch(GenomicSearch):
    def __init__(self):
        self._features = collections.defaultdict(
            lambda: collections.defaultdict(FeatureSet)
        )
        self._starts = collections.defaultdict(list)
        self._locked: bool = False

    def add_feature(self, location: Location, feature: Any):
        self._locked = False

        if (
            location.chr not in self.features
            or location.start not in self.features[location.chr]
        ):
            self.features[location.chr][location.start] = FeatureSet(location.start)

        if (
            location.chr not in self.features
            or location.end not in self.features[location.chr]
        ):
            self.features[location.chr][location.end] = FeatureSet(location.end)

        self._features[location.chr][location.start].add(feature)
        self._features[location.chr][location.end].add(feature)

    def _sort_starts(self):
        """
        Sort the starts for binary searching
        """

        if self._locked:
            return

        # We need the locations sorted in order for a binary search
        for chr in self._features:
            self._starts[chr] = sorted(self.features[chr])

        self._locked = True

    def get_features(self, location: Location) -> list[FeatureSet]:
        """
        Return a list of features closest to a location which then
        be tested for overlap etc
        """

        ret = []

        if location.chr not in self._features:
            return ret

        self._sort_starts()

        starts = self._starts[location.chr]

        s = self._get_start_index(starts, location.start, True)
        e = self._get_start_index(starts, location.end, False)

        chr_features = self._features[location.chr]

        for i in range(s, e + 1):
            ret.append(chr_features[starts[i]])

        return ret

    def _get_start_index(self, indices: list[int], start: int, startMode: bool):
        if len(indices) < 2:
            return 0

        s = 0

        if start <= indices[s]:
            return s

        e = len(indices) - 1

        if start >= indices[e]:
            return e

        while e - s > 1:
            im = int((s + e) / 2)  # s + (e - s) // 2

            pm = indices[im]

            if pm > start:
                e = im  # - 1
            elif pm < start:
                s = im  # + 1
            else:
                return im

        # If we made it this far, we have narrowed down the search to
        # two position so return based on trying to cover the region
        # of interest. We can narrow down the actual closest bin later.

        if startMode:
            # return the closest before the start
            return s
        else:
            # return the closest after the end
            return e


class BlockSearch(GenomicSearch):
    def __init__(self, block_size: int = 10000):
        self._block_size = block_size
        self._max_bin = 0
        self._features = collections.defaultdict(
            lambda: collections.defaultdict(FeatureSet)
        )

    def add_feature(self, location: Location, feature: Location):
        s = int(location.start / self._block_size)
        e = int(location.end / self._block_size)

        for b in range(s, e + 1):
            if (
                location.chr not in self._features
                or b not in self._features[location.chr]
            ):
                self._features[location.chr][b] = FeatureSet(location.start)

            self._features[location.chr][b].add(feature)
            self._max_bin = max(self._max_bin, b)

    # def get_closest_features(self, location: Location) -> FeatureSet:
    #     featureSets = self.get_features(location)

    #     min_d = sys.maxint
    #     ret = None

    #     for featureSet in featureSets:
    #         d = abs(mid(location) - featureSet.start)

    #         if d < min_d:
    #             ret = featureSet
    #             min_d = d

    #     return ret

    def get_closest_features(self, location: Location) -> FeatureSet:
        if location.chr not in self._features:
            return None

        s = int(location.start / self._block_size)
        e = int(location.end / self._block_size)

        chr_features = self._features[location.chr]

        # extend s and and e to closest bins with actual features since
        # we want the closest

        while s > 0:
            if s in chr_features:
                break
            s -= 1

        while e < self._max_bin:
            if e in chr_features:
                break
            e += 1

        featureSets: list[FeatureSet] = []

        for b in range(s, e + 1):
            if b in chr_features:
                featureSets.append(chr_features[b])

        min_d = sys.maxsize
        ret = None

        for featureSet in featureSets:
            d = abs(mid(location) - featureSet.start)

            if d < min_d:
                ret = featureSet
                min_d = d

        return ret

    def get_features(self, location: Location) -> list[FeatureSet]:
        """
        Return a list of features closest to a location which then
        be tested for overlap etc

        Args:
                location (Location): location to annotate

        Returns:
                list[FeatureSet]: set of features for testing for strict overlap
        """

        ret: list[FeatureSet] = []

        if location.chr not in self._features:
            return ret

        s = int(location.start / self._block_size)
        e = int(location.end / self._block_size)

        chr_features = self._features[location.chr]

        # extend s and

        for b in range(s, e + 1):
            if b in chr_features:
                ret.append(chr_features[b])

        return ret


def is_location(location: str):
    """
    Returns true if string looks like a location

    Args:
            location (str): a location string, e.g "chr1:1-2"

    Returns:
            bool: True if string looks like a location
    """
    return re.match(r"(chr.+):(\d+)-(\d+)", location)


def is_chr(location: str) -> bool:
    """
    Returns true if location is a chr

    Args:
            location (str): a location string, e.g "chr1:1-2"

    Returns:
            bool: True if contains "chr"
    """
    return re.match(r"(chr.+)", location)


def parse_location(location: str, padding5p: int = 0, padding3p: int = 0) -> Location:
    """
    Parses a location string, e.g. 'chr1:1-100' into a location object.

    Args:
            location (str): A location string
            padding5p (int, optional): Add padding to the 5p end. Defaults to 0.
            padding3p (int, optional): Add padding to the 3p end. Defaults to 0.

    Returns:
            Location: A location object representing the location string.
    """

    matcher = re.match(
        r".*(chr.+)-(\d+)-(\d+).*", location.replace(",", "").replace(":", "-")
    )

    chr = matcher.group(1)
    start = int(matcher.group(2))
    end = int(matcher.group(3))

    return Location(chr, start - padding5p, end + padding3p)


def parse_location_cols(tokens: list[str], offset: int = 0) -> Location:
    return Location(tokens[offset], int(tokens[offset + 1]), int(tokens[offset + 2]))


def pad_location(location: Location, padding5p: int = 0, padding3p: int = 0):
    return Location(location.chr, location.start - padding5p, location.end + padding3p)


def max_region(location1: Location, location2: Location) -> Union[None, Location]:
    if location1.chr != location2.chr:
        return None

    min_start = min(location1.start, location2.start)
    # max_end = max(location1.end, location2.end)
    # min_start = min(location1.start, location2.start)
    max_end = max(location1.end, location2.end)

    return Location(location1.chr, min_start, max_end)


def overlap_locations(
    location1: Location, location2: Location
) -> Union[None, Location]:
    """Overlap two locations and return a location representing the overlap

    Args:
            location1 (Location): _description_
            location2 (Location): _description_

    Returns:
            Union[None, Location]: _description_
    """
    if location1.chr != location2.chr:
        return None

    max_start = max(location1.start, location2.start)
    # max_end = max(location1.end, location2.end)
    # min_start = min(location1.start, location2.start)
    min_end = min(location1.end, location2.end)

    if min_end < max_start:
        return None

    # if location1.end < location2.start or \
    #         location2.end < location1.start or \
    #         location1.start > location2.end or \
    #         location2.start > location1.end:
    #     return None

    # start = max(location1.start, location2.start)
    # end = min(location1.end, location2.end)

    return Location(location1.chr, max_start, min_end)


def max_region_from_locs(locs: list[Location]) -> Union[Location, None]:
    """
    Given a list of locations, return the maximum overlap
    """

    # locs = list(sorted(locs))

    start = min([loc.start for loc in locs])
    end = max([loc.end for loc in locs])

    return Location(locs[0].chr, start, end)


def overlap_fraction(location1: Location, location2: Location) -> float:
    """Fraction that location1 that overlaps with location2 (denominator)

    Args:
            location1 (Location): _description_
            location2 (Location): _description_

    Returns:
            Union[None, Location]: _description_
    """
    if location1 is None or location2 is None:
        return 0

    if location1.chr != location2.chr:
        return 0

    max_start = max(location1.start, location2.start)
    # max_end = max(location1.end, location2.end)
    # min_start = min(location1.start, location2.start)
    min_end = min(location1.end, location2.end)

    if min_end < max_start:
        return 0

    return (min_end - max_start + 1) / location2.length


def is_overlapping(location1: Location, location2: Location) -> bool:
    if location1.chr != location2.chr:
        return False

    max_start = max(location1.start, location2.start)
    # max_end = max(location1.end, location2.end)
    # min_start = min(location1.start, location2.start)

    min_end = min(location1.end, location2.end)

    return max_start <= min_end  # return max_start >= min_start and max_start < min_end


def get_closest_tss(tss):
    closest = NA

    min_d = 1000000
    min_abs_d = 1000000

    for p in tss:
        if p == "n/a" or p == "NA":
            continue

        d = int(p)
        abs_d = abs(d)

        if abs_d < min_abs_d:
            min_d = d
            min_abs_d = abs_d

    if min_d != 1000000:
        closest = str(min_d)

    return closest


def load_locations(file, center=False, padding5p=0, padding3p=0):
    """
    Load a set of locations
    """

    f = open(file, "r")

    header = f.readline().strip()

    ret = []

    for line in f:
        l = parse_location(line.strip())

        if center:
            l = center_location(l)

        l = pad_location(l, padding5p, padding3p)

        ret.append(l)

    f.close()

    return ret, header


class SearchGenomicFeatures(GappedSearch):
    """
    Instance of gapped search for quickly identifying if a region
    is close to centromere
    """

    def __init__(self, file: str):
        super().__init__()

        f = open(file, "r")

        f.readline()

        for line in f:
            line = line.strip()

            if len(line) == 0:
                continue

            location: Location = parse_location(line)

            # we store locations for determining if they overlap or not
            self.add_feature(location, location)

        f.close()

        self.lock()


class GenomicFeaturesOverlap:
    """
    Uses a gapped search to determine by how much a location overlaps a feature
    """

    def __init__(self, search: BlockSearch):
        self._search = search

    def get_features(self, location: Location):
        return self._search.get_features(location)

    def get_max_overlap(self, location: Location):
        features = self.get_features(location)

        max_overlap_width = -1
        max_overlap = None

        for feature in features:
            for l in feature.values:
                overlap = overlap_locations(location, l)

                if overlap is not None:
                    if overlap.length > max_overlap_width:
                        max_overlap_width = overlap.length
                        max_overlap = overlap

        return max_overlap


class Annotation(ABC):
    """
    Modules for annotating a peak by appending columns to a row of
    peak meta data.
    """

    def __init__(
        self,
        name: str,
        species: species.Species = species.Species.HUMAN,
        version: str = "1.0.0",
    ):
        self._metadata = {
            "name": name,
            "species": species.name.lower(),
            "version": version,
        }

    @abstractmethod
    def get_names(self) -> list[str]:
        """
        Returns the heading names to add.

        Returns:
                list[str]: headers.
        """
        ...

    # def annotate(self, location: Location) -> list[Union[str, int, float]]:
    #     return []

    @abstractmethod
    def update_row(self, location: Location, row_map: Mapping[str, Any]) -> list[Any]:
        """
        Given existing row annotation, return the new columns to add. Columns should
        match the headers returned in get_names().

        Args:
                location (Location): Location to annotated.
                row_map (Mapping[str, Any]):    Annotations already assigned using
                                                                                header name as key. This is to
                                                                                allow annotation to be based on
                                                                                previous modules.

        Returns:
                list[Any]: list of columns.
        """
        ...

    @property
    def metadata(self) -> dict[str, Any]:
        """
        Returns info about the module and how it was run for logging purposes to
        keep a record of how a peak table was annotated.

        Returns:
                dict[str, Any]: dictionary of keys describing module.
        """

        return self._metadata


class ClassifyRegion(ABC):
    @abstractmethod
    def get_classification(self, location: Location) -> str:
        """
        Returns a class label for a region, such as identifying a telomere.

        Args:
                location (Location): Location to classify.

        Returns:
                str: a class label.
        """
        ...

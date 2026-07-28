"""
This module provides utilities for processing databases, parsing sequences, checking different file or object formats and logging.
"""  # noqa: EXE002

from ._checkers_and_helpers import (
    CustomYamlDumper,
    cast_to_int,
    cast_to_list,
    cast_to_list_of_lists,
    cast_to_string,
    check_if_dna_sequence,
    check_if_key_exists,
    check_tsv_format,
    generate_unique_filename,
)
from ._database_processor import (
    check_if_key_in_database,
    check_if_region_in_database,
    collapse_properties_for_duplicated_sequences,
    flatten_property_list,
    format_oligo_properties,
    merge_databases,
)
from ._logging import configure_root_logger, logger
from ._sequence_parser import FastaParser, GffParser, VCFParser
from ._sequence_processor import (
    append_nucleotide_to_sequences,
    count_kmer_abundance,
    get_complement_regions,
    get_intersection,
    get_sequence_from_annotation,
    remove_index_files,
)

__all__ = [
    "CustomYamlDumper",
    "FastaParser",
    "GffParser",
    "VCFParser",
    "append_nucleotide_to_sequences",
    "cast_to_int",
    "cast_to_list",
    "cast_to_list_of_lists",
    "cast_to_string",
    "check_if_dna_sequence",
    "check_if_key_exists",
    "check_if_key_in_database",
    "check_if_region_in_database",
    "check_tsv_format",
    "collapse_properties_for_duplicated_sequences",
    "configure_root_logger",
    "count_kmer_abundance",
    "flatten_property_list",
    "format_oligo_properties",
    "generate_unique_filename",
    "get_complement_regions",
    "get_intersection",
    "get_sequence_from_annotation",
    "logger",
    "merge_databases",
    "remove_index_files",
]

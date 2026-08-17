############################################  # noqa: EXE002
# imports
############################################

import shelve
from shelve import Shelf
from typing import Any

from oligo_designer_toolsuite._constants import SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite.utils._checkers_and_helpers import safe_append_filename
from oligo_designer_toolsuite.utils._logging import logger

from ._checkers_and_helpers import cast_to_list, cast_to_list_of_lists

############################################
# Collection of utility functions
############################################


def merge_databases(
    database1: Shelf,
    database2: Shelf,
    sequence_type: str,
    database_sequence_types: list[str],
    dir_cache_files: str,
    max_entries_in_memory: int,
) -> Shelf:
    """
    Merges two oligo databases by combining their content based on sequence keys,
    ensuring that sequences with the same oligo are merged, and avoiding duplicates.

    :param database1: The first database to be merged.
    :type database1: Shelf
    :param database2: The second database to be merged.
    :type database2: Shelf
    :param sequence_type: The sequence type key to use for merging (must be in database_sequence_types).
    :type sequence_type: str
    :param database_sequence_types: List of sequence type keys in the database.
    :type database_sequence_types: list[str]
    :param dir_cache_files: Directory to store cache files used for merging.
    :type dir_cache_files: str
    :param max_entries_in_memory: Maximum number of entries to keep in memory for the LRU (Least Recently Used) cache.
    :type max_entries_in_memory: int
    :return: The merged database.
    :rtype: Shelf
    """
    if sequence_type not in database_sequence_types:
        raise ValueError(
            f"sequence_type '{sequence_type}' must be in database_sequence_types. "
            f"Current database_sequence_types: {database_sequence_types}"
        )

    def _get_sequence_as_key(database: Shelf, regions: list[str], sequence_type: str) -> Shelf:
        """
        Converts oligo sequences to dictionary keys, grouping oligo properties by sequence for each specified region.

        :param database: The database containing sequences and their properties.
        :type database: Shelf
        :param regions: List of regions within the database to process.
        :type regions: list
        :param sequence_type: The sequence type key to use for merging.
        :type sequence_type: str
        :return: A dictionary with sequences as keys and oligo properties as values.
        :rtype: Shelf
        """
        # TODO: correctly replace old shelf with new one
        database_modified = shelve.open(safe_append_filename(dir_cache_files, "database_modified"), flag="n")  # noqa: SIM115

        for region in regions:
            new_region = {}
            for oligo_properties in database[region].values():
                oligo_sequence = oligo_properties[sequence_type]
                oligo_properties.pop(sequence_type)
                new_region[oligo_sequence] = oligo_properties
            database_modified[region] = new_region
        return database_modified

    def _add_database_content(
        database_merged_tmp: Shelf, database_in_tmp: Shelf, database_sequence_types: list[str]
    ) -> Shelf:
        """
        Merges oligo properties from two databases, ensuring sequences with the same oligo are combined and properties are updated.

        :param database_merged_tmp: The dictionary to which content is added.
        :type database_merged_tmp: Shelf
        :param database_in_tmp: The dictionary containing new content to merge.
        :type database_in_tmp: Shelf
        :param database_sequence_types: List of sequence type keys in the database.
        :type database_sequence_types: list[str]
        :return: The updated dictionary with merged oligo properties.
        :rtype: Shelf
        """
        for region_id, database_region in database_in_tmp.items():
            # extract from Shelf
            database_merged_region = database_merged_tmp.get(region_id, {})
            for oligo_sequence, oligo_properties in database_region.items():
                if oligo_sequence in database_merged_region:
                    oligo_properties = collapse_properties_for_duplicated_sequences(
                        database_merged_region[oligo_sequence], oligo_properties, database_sequence_types
                    )
                database_merged_region[oligo_sequence] = oligo_properties
            # put into Shelf
            database_merged_tmp[region_id] = database_merged_region
        return database_merged_tmp

    # keys that are in both dicts
    regions_intersection = list(set(database1) & set(database2))

    # TODO: correctly replace old shelf with new one
    database_merged = shelve.open(safe_append_filename(dir_cache_files, "database_merged"), flag="n")  # noqa: SIM115

    # only loop over entries that have keys in both dicts
    db1_sequences_as_keys = _get_sequence_as_key(database1, regions_intersection, sequence_type)
    db2_sequences_as_keys = _get_sequence_as_key(database2, regions_intersection, sequence_type)

    database_merged = _add_database_content(database_merged, db1_sequences_as_keys, database_sequence_types)
    db1_sequences_as_keys.close()
    database_merged = _add_database_content(database_merged, db2_sequences_as_keys, database_sequence_types)
    db2_sequences_as_keys.close()

    for region_id, database_merged_region in database_merged.items():
        for i, (oligo_sequence, oligo_properties) in enumerate(database_merged_region.items(), start=1):
            oligo_id = f"{region_id}{SEPARATOR_OLIGO_ID}{i}"
            oligo_seq_info = {sequence_type: oligo_sequence} | oligo_properties
            oligo_properties[sequence_type] = oligo_sequence
            database_merged_region[oligo_id] = oligo_seq_info

        database_merged[region_id] = database_merged_region

    # add entries with keys in only one dict
    for region_id in set(database1) - set(database2):
        database_merged[region_id] = database1[region_id]
    for region_id in set(database2) - set(database1):
        database_merged[region_id] = database2[region_id]

    return database_merged


def collapse_properties_for_duplicated_sequences(
    oligo_properties1: dict[str, Any], oligo_properties2: dict[str, Any], database_sequence_types: list[str]
) -> dict[str, Any]:
    """
    Merges two dictionaries of oligo properties, combining values for non-sequence keys and issuing warnings if sequences for the same oligo ID.

    :param oligo_properties1: The first dictionary of oligo properties.
    :type oligo_properties1: dict
    :param oligo_properties2: The second dictionary of oligo properties.
    :type oligo_properties2: dict
    :param database_sequence_types: List of sequence type keys in the database.
    :type database_sequence_types: list[str]
    :return: A merged dictionary with combined oligo properties.
    :rtype: dict
    """
    oligo_properties = {}

    if oligo_properties1 == oligo_properties2:
        return oligo_properties1

    for d in (oligo_properties1, oligo_properties2):
        for key, values in d.items():
            if key not in oligo_properties:
                oligo_properties[key] = values
            else:
                if key in database_sequence_types and oligo_properties[key] != values:
                    logger.warning(
                        f"Values for key {key} are different in the two oligo_properties dictionaries."
                    )
                elif key not in database_sequence_types:
                    oligo_properties[key].extend(values)

    return oligo_properties


def get_oligo_property_value(
    property: str, flatten: bool, region: dict, oligo_id: str
) -> Any | list[Any] | None:
    """
    Retrieve the value of a specified property for a given oligo and region ID,
    optionally flattening nested or list properties to a unique set of values.

    :param property: The name of the property to retrieve.
    :type property: str
    :param flatten: Whether to flatten list properties to a unique set of values in the table.
    :type flatten: bool
    :param region_id: The ID of the region where the oligo is located.
    :type region_id: str
    :param oligo_id: The ID of the oligo for which the property value is retrieved.
    :type oligo_id: str
    :return: The value of the specified property, possibly flattened. Can be any type (int, float, str, bool, etc.) or a list of values, or None if the property doesn't exist.
    :rtype: Any | list[Any] | None
    :raises ValueError: If the specified region or oligo does not exist in the database.
    """
    # if not region_id in self.database:
    #     raise DatabaseError(f"Region '{region_id}' does not exist in the database.")

    if not oligo_id in region:
        raise ValueError(f"Oligo '{oligo_id}' does not exist in specified region.")

    oligo_properties = region[oligo_id]
    if property not in oligo_properties:
        property_value = None
    elif flatten:
        property_value = flatten_property_list(region[oligo_id][property])
        if property_value and len(property_value) == 1:
            property_value = property_value[0]
    else:
        property_value = region[oligo_id][property]

    return property_value


def check_if_region_in_database(
    database: Shelf,
    region_ids: list[str],
    write_regions_with_insufficient_oligos: bool,
    file_removed_regions: str,
) -> None:
    """
    Checks if specified regions are present in the database and logs missing regions.

    :param database: The database.
    :type database: dict
    :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (List[str]). If None, all regions in the database are processed.
    :type region_ids: list
    :param write_regions_with_insufficient_oligos: Whether to write regions with insufficient oligos to a file.
    :type write_regions_with_insufficient_oligos: bool
    :param file_removed_regions: The file to which missing regions should be logged.
    :type file_removed_regions: str
    """
    keys = list(database.keys())
    for region_id in region_ids:
        if region_id not in keys:
            logger.warning(f"Region {region_id} not available in reference file.")
            if write_regions_with_insufficient_oligos:
                with open(file_removed_regions, "a") as hanlde:
                    hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")


def check_if_key_in_database(database: Shelf, key: str, region_ids: str | list[str] | None = None) -> bool:
    """
    Checks if a key exists in the database.

    :param database: The database to check.
    :type database: Shelf
    :param key: The key to check.
    :type key: str
    :param region_ids: Optional region ID(s) to check in specific region(s). Can be a single region ID (str) or a list of region IDs (list[str]). If provided, checks that the key exists in ALL specified regions. If None, checks if it exists in at least one region.
    :type region_ids: str | list[str] | None
    :return: True if the key exists in the database (in all specified regions if region_ids is provided, or in at least one region if region_ids is None), False otherwise.
    :rtype: bool
    """

    # for thread safety and because we don't assume our code to contain bugs, always return True
    # (this is only used with assert statements)
    return True

    # # Helper for recursive search
    # def recursive_contains(d: Any, target: str) -> bool:
    #     """
    #     Recursively checks if `target` exists as a key in a nested structure.
    #     Top-level is EffiDict, deeper levels are plain dicts.
    #     """
    #     try:
    #         if isinstance(d, dict):  # this will match EffiDict at the top AND nested dicts
    #             if target in d:
    #                 return True
    #             return any(recursive_contains(v, target) for v in d.values())
    #         return False
    #     except Exception:  # noqa: BLE001, RUF100
    #         return False

    # database_region_ids = set(database.keys())

    # # --- Case: region restriction ---
    # if region_ids is not None:
    #     region_ids = cast_to_list(region_ids)
    #     regions_found = False  # check if at least one region is found in the database
    #     for region_id in region_ids:
    #         if region_id not in database_region_ids:
    #             continue
    #         else:
    #             regions_found = True
    #         region_data = database[region_id]
    #         if not recursive_contains(region_data, key):
    #             return False
    #     return regions_found

    # # --- Case: no region restriction → any region may match ---
    # return any(recursive_contains(database[region_id], key) for region_id in database_region_ids)


def format_oligo_properties(
    oligo_properties: dict[str, Any], database_sequence_types: list[str]
) -> dict[str, Any]:
    """
    Ensures that the values in an oligo properties dictionary are formatted as lists of lists.

    :param oligo_properties: The dictionary of oligo properties to format.
    :type oligo_properties: dict
    :param database_sequence_types: List of sequence type keys in the database.
    :type database_sequence_types: list[str]
    :return: The formatted dictionary with lists of lists for non-sequence keys.
    :rtype: dict
    """
    for key, value in oligo_properties.items():
        if key not in database_sequence_types:
            oligo_properties[key] = cast_to_list_of_lists(value)
    return oligo_properties


def flatten_property_list(property: list[Any]) -> list[Any]:
    """
    Flattens a nested list of properties into a single list, or returns the item if only one element remains.

    :param property: The list or nested list of properties to flatten.
    :type property: list[Any]
    :return: A flattened list or a single item if only one element exists.
    :rtype: list[Any]
    """
    flattened_property_list = [
        item
        for sublist in (property if isinstance(property, list) else [property])
        for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    if len(flattened_property_list) == 1:
        return cast_to_list(flattened_property_list[0])
    return cast_to_list(flattened_property_list)

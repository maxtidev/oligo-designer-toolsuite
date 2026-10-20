############################################  # noqa: EXE002
# imports
############################################

from typing import Any

import numpy as np
import pandas as pd

from oligo_designer_toolsuite.utils._logging import logger

from ._checkers_and_helpers import cast_to_list

############################################
# Collection of utility functions
############################################


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
    property: str, flatten: bool, region: pd.DataFrame, oligo_id: str
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
    if property not in region:
        return None
    try:
        property_value = region.at[oligo_id, property]
    except KeyError:
        raise ValueError(f"Oligo '{oligo_id}' does not exist in specified region.")

    if property_value is pd.NA:
        return None

    if flatten:
        property_value = flatten_property_list(property_value)

        if property_value and len(property_value) == 1:
            property_value = property_value[0]

        # Due to some functions not handling lists of None values properly, merge lists of None into just None
        # also see exception 2 in load_database_from_parquet
        if isinstance(property_value, list) and all(v is None for v in property_value):
            property_value = None

    return property_value


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
            # We assume these values to be lists and cast them to lists of lists
            oligo_properties[key] = [value]
    return oligo_properties


def flatten_property_list(property: Any) -> list[Any]:
    """
    Flattens a nested list of properties into a single list, or returns the item if only one element remains.

    :param property: The list or nested list of properties to flatten.
    :type property: list[Any]
    :return: A flattened list or a single item if only one element exists.
    :rtype: list[Any]
    """
    flattened_property_list = [
        item
        for sublist in (property if isinstance(property, list | np.ndarray) else [property])
        for item in (sublist if isinstance(sublist, list | np.ndarray) else [sublist])
    ]
    if len(flattened_property_list) == 1:  # NOTE: questionable
        return cast_to_list(flattened_property_list[0])
    return cast_to_list(flattened_property_list)

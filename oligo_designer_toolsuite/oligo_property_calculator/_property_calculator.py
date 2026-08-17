############################################
# imports
############################################

from typing import Any

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import BaseProperty
from oligo_designer_toolsuite.utils import check_if_key_in_database

############################################
# Property Calculator Class
############################################


class PropertyCalculator:
    """
    A class for applying multiple property calculators to oligonucleotides in an OligoDatabase.

    The `PropertyCalculator` class allows you to apply a list of property calculators (subclasses of `BaseProperty`) to an OligoDatabase.
    The properties are calculated in parallel across all regions of the database, and the calculated values are stored as properties in the database.

    :param properties: A list of property calculators to apply to oligonucleotides.
    :type properties: list[BaseProperty]
    """

    def __init__(self, properties: list[BaseProperty]) -> None:
        """Constructor for the PropertyCalculator class."""
        self.properties = properties

    def apply(self, oligo_database: OligoDatabase, sequence_type: str, n_jobs: int = 1) -> OligoDatabase:
        """
        Apply the property calculators to all oligonucleotides in the OligoDatabase and update
        the database with the calculated property values.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param n_jobs: Number of parallel jobs to use for processing. Defaults to 1.
        :type n_jobs: int
        :return: The updated OligoDatabase with the calculated properties.
        :rtype: OligoDatabase
        """
        assert check_if_key_in_database(oligo_database.database, sequence_type), (
            f"Sequence type '{sequence_type}' not found in database."
        )

        oligo_database.map_regions(
            self._calculate_region,
            args=(sequence_type,),
            description="Property Calculator",
            n_jobs=n_jobs,
        )

        return oligo_database

    def _calculate_region(self, oligo_database: OligoDatabase, region_id: str, sequence_type: str) -> None:
        """
        Calculate properties for all oligonucleotides in a specific region of the OligoDatabase.

        This method iterates through the oligonucleotides in a given region of the database,
        applying all property calculators to each oligo and updating the database with the calculated values.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        """
        database_region = oligo_database.load_region(region_id)
        new_oligo_property: dict[str, dict[str, Any]] = {}

        for oligo_id in database_region.keys():  # noqa: SIM118
            # Calculate all properties for this oligo
            for property_calc in self.properties:
                property_result = property_calc.apply(
                    region=database_region,
                    oligo_id=oligo_id,
                    sequence_type=sequence_type,
                )
                # Merge results into the property dictionary
                if oligo_id not in new_oligo_property:
                    new_oligo_property[oligo_id] = {}
                new_oligo_property[oligo_id].update(property_result)

        # Update only this region (avoids O(regions × total_oligos) full-database scan per worker)  # noqa: RUF003
        oligo_database.update_oligo_properties(new_oligo_property, region_ids=region_id)

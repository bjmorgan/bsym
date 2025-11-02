from __future__ import annotations

import numpy as np
from bsym.symmetry_group import SymmetryGroup
from bsym.configuration_space import ConfigurationSpace


class CoordinateConfigSpace(ConfigurationSpace):
    """
    A :any:`CoordinateConfigSpace` object is a :any:`ConfigurationSpace` that has an associated
    set of coordinates. Each vector in the configuration vector space has a corresponding coordinate.
    """
    
    def __init__(
        self,
        coordinates: np.ndarray,
        symmetry_group: SymmetryGroup | None = None,
        objects: np.ndarray | list[int] | None = None
    ) -> None:
        """
        Create a :any:`CoordinateConfigSpace` object.
        
        Args:
            coordinates: The set of coordinates that describe the vector space of this 
                configuration space.
            symmetry_group: The set of symmetry operations describing the symmetries of 
                this configuration space.
            objects: Optional array or list of objects to represent the coordinates. If None, 
                creates objects as np.arange(len(coordinates)) + 1.
        
        Returns:
            None
        """
        if objects is None:
            # Create a set of objects to represent the coordinates.
            objects = np.arange(len(coordinates)) + 1
        # Convert to list for parent class
        objects_list = objects.tolist() if isinstance(objects, np.ndarray) else list(objects)
        super().__init__(objects_list, symmetry_group)
        self.coordinates = coordinates
    
    def unique_coordinates(
        self,
        site_distribution: dict[int, int],
        verbose: bool = False
    ) -> list[dict[int, np.ndarray]]:
        """
        Find the symmetry inequivalent coordinates for a given site occupation.
        
        Args:
            site_distribution: A dictionary that defines the number of each object 
                to be arranged in this system.
                e.g. for a structure with four sites, with two occupied (denoted `1`)
                and two unoccupied (denoted `0`)::
                    {1: 2, 0: 2}
            verbose: Print verbose output.
        
        Returns:
            A list of dicts. Each dict describes the set of coordinates for each site type.
        """
        unique_configs = self.unique_configurations(site_distribution, verbose=verbose)
        unique_coordinates = [u.map_objects(self.coordinates) for u in unique_configs]
        return unique_coordinates
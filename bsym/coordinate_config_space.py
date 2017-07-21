from bsym import ConfigurationSpace
import numpy as np

class CoordinateConfigSpace( ConfigurationSpace ):
    """
    A :any:`CoordinateConfigSpace` object is a :any`ConfigurationSpace` that has an associated
    set of coordinates. Each vector in the configuration vector space has a corresponding coordinate.
    """

    def __init__( self, coordinates, symmetry_group=None, objects=None ):
        """
        Create a :any:`CoordinateConfigSpace` object.

        Args:
            coordinates (np.array): The set of coordinates that describe the vector space of this configuration space.
            symmetry_group (:any:`SymmetryGroup`): The set of symmetry operations describing the symmetries of this configuration space.

        Returns:
            None
        """
        if objects is None:
            # Create a set of objects to represent the coordinates.
            objects = np.arange( len( coordinates ) ) + 1
        super().__init__( objects, symmetry_group )
        self.coordinates = coordinates

    def unique_coordinates( self, site_distribution, verbose=False ):
        """
        Find the symmetry inequivalent coordinates for a given site occupation.

        Args:
            site_distribution (dict): A dictionary that defines the number of each object 
                                      to be arranged in this system.

                                      e.g. for a structure with four sites, with two occupied (denoted `1`)
                                      and two unoccupied (denoted `0`)::

                                          { 1 : 2, 2 : 1 }
            verbose (opt:default=False): Print verbose output.

        Returns:
            unique_coordinates (list[dict]): A list of dicts. Each dict describes the set of coordinates for each site type.
        """
        unique_configs = self.unique_configurations( site_distribution, verbose=verbose )
        unique_coordinates = [ u.map_objects( self.coordinates ) for u in unique_configs ]
        return unique_coordinates

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.util.coord_utils import coord_list_mapping_pbc
from pymatgen import Lattice, Structure

from bsym.spacegroup import SpaceGroup
from bsym import symmetry_operation as so

def unique_symmetry_operations_as_vectors_from_structure( structure ):
    """
    Uses pymatgen symmetry analysis to find the minimum complete set of symmetry operations for the space group of a structure.

    Args:
        structure (pymatgen Structure): structure to be analysed.

    Returns:
        a list of lists, containing the symmetry operations as vector mappings.
    """
    symmetry_analyzer = SpacegroupAnalyzer( structure )
    symmetry_operations = symmetry_analyzer.get_symmetry_operations()

    mappings = []
    for symmop in symmetry_operations:
        new_structure = Structure( structure.lattice, structure.species, symmop.operate_multi( structure.cart_coords ) )
        new_mapping = [ x+1 for x in list( coord_list_mapping_pbc( new_structure.cart_coords, structure.cart_coords ) ) ]
        if new_mapping not in mappings:
            mappings.append( new_mapping )

    return mappings

def spacegroup_from_structure( structure ):
    """
    Generates a SpaceGroup object from a pymatgen Structure 
    Args:
        structure (pymatgen Structure): structure to be used to define the SpaceGroup.

    Returns:
        a new SpaceGroup instance
    """
    mappings = unique_symmetry_operations_as_vectors_from_structure( structure )
    symmetry_operations = [ so.SymmetryOperation.from_vector( m ) for m in mappings ]
    return SpaceGroup( symmetry_operations = symmetry_operations )

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.util.coord_utils import coord_list_mapping_pbc
from pymatgen import Lattice, Structure

from bsym.spacegroup import SpaceGroup
from bsym import symmetry_operation as so

def unique_symmetry_operations_as_vectors_from_structure( structure, subset = None ):
    """
    Uses pymatgen symmetry analysis to find the minimum complete set of symmetry operations for the space group of a structure.

    Args:
        structure (pymatgen Structure): structure to be analysed.
        subset    (Optional [list]):    list of atom indices to be used for generating the symmetry operations

    Returns:
        a list of lists, containing the symmetry operations as vector mappings.
    """
    symmetry_analyzer = SpacegroupAnalyzer( structure )

    print( "The spacegroup for structure is {}".format(symmetry_analyzer.get_spacegroup_symbol()) )

    symmetry_operations = symmetry_analyzer.get_symmetry_operations()

    mappings = []
    if subset:
        species_subset = [ spec for i,spec in enumerate(structure.species) if i in subset]
        frac_coords_subset = [ coord for i, coord in enumerate( structure.frac_coords ) if i in subset ]
        mapping_structure = Structure( structure.lattice, species_subset, frac_coords_subset ) 
    else:
        mapping_structure = structure

    for symmop in symmetry_operations:
        new_structure = Structure( mapping_structure.lattice, mapping_structure.species, symmop.operate_multi( mapping_structure.frac_coords ) )
        new_mapping = [ x+1 for x in list( coord_list_mapping_pbc( new_structure.frac_coords, mapping_structure.frac_coords ) ) ]
        if new_mapping not in mappings:
            mappings.append( new_mapping )

    return mappings

def spacegroup_from_structure( structure, subset = None ):
    """
    Generates a SpaceGroup object from a pymatgen Structure 
    Args:
        structure (pymatgen Structure): structure to be used to define the SpaceGroup.
        subset    (Optional [list]):    list of atom indices to be used for generating the symmetry operations
    Returns:
        a new SpaceGroup instance
    """
    mappings = unique_symmetry_operations_as_vectors_from_structure( structure, subset )
    symmetry_operations = [ so.SymmetryOperation.from_vector( m ) for m in mappings ]
    return SpaceGroup( symmetry_operations = symmetry_operations )

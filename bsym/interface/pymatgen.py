from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.util.coord_utils import coord_list_mapping_pbc
from pymatgen import Lattice, Structure

from bsym import SpaceGroup
from bsym import SymmetryOperation

def unique_symmetry_operations_as_vectors_from_structure( structure, verbose=True, subset=None ):
    """
    Uses `pymatgen`_ symmetry analysis to find the minimum complete set of symmetry operations for the space group of a structure.

    Args:
        structure (pymatgen ``Structure``): structure to be analysed.
        subset    (Optional [list]):        list of atom indices to be used for generating the symmetry operations.

    Returns:
        (list[list]): a list of lists, containing the symmetry operations as vector mappings.

    .. _pymatgen:
        http://pymatgen.org

    """
    symmetry_analyzer = SpacegroupAnalyzer( structure )
    if verbose:
        print( "The spacegroup for this structure is {}".format( symmetry_analyzer.get_space_group_symbol()) )
    symmetry_operations = symmetry_analyzer.get_symmetry_operations()
    mappings = []
    if subset:
        species_subset = [ spec for i,spec in enumerate( structure.species ) if i in subset]
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
    Generates a ``SpaceGroup`` object from a `pymatgen` ``Structure``. 

    Args:
        structure (pymatgen ``Structure``): structure to be used to define the :any:`SpaceGroup`.
        subset    (Optional [list]):        list of atom indices to be used for generating the symmetry operations.

    Returns:
        a new :any:`SpaceGroup` instance
    """
    mappings = unique_symmetry_operations_as_vectors_from_structure( structure, subset )
    symmetry_operations = [ SymmetryOperation.from_vector( m ) for m in mappings ]
    return SpaceGroup( symmetry_operations=symmetry_operations )

def poscar_from_sitelist( configs, labels, sitelists, structure, subset=None ):
    """
    Uses `pymatgen` ``Structure.to()`` method to generate `VASP` `POSCAR` files for a set of 
    configurations within a parent structure.

    Args:
        configs   (list [:any:`Configuration`]): list of :any:`configurations`.
        labels    (list [int]):                  labels defining order of sites in output.
        sitelists (list [:any;`SiteList`]):      list of :any:`SiteList` objects.
        structure (pymatgen ``Structure``):      parent structure.
        subset    (Optional [list]):             list of atom indices to include from the parent structure. 

    Returns:
        None
    """
    if subset:
        species_clean = [ spec for i,spec in enumerate( structure.species ) if i not in subset ]
        species_config = [ spec for i,spec in enumerate( structure.species ) if i in subset ]
        frac_coords_clean = [ coord for i, coord in enumerate( structure.frac_coords ) if i not in subset ]
        clean_structure = Structure( structure.lattice, species_clean, frac_coords_clean )
    else:
        clean_structure = Structure( structure.lattice, [], [] )
        species_config = structure.species
    for idx, config in enumerate( configs, start=1 ):
       structure_config = clean_structure.copy()
       for label in labels:
           for pos in config.position( label ):
               for sitelist in sitelists:
                   structure_config.append( species_config[ pos ], sitelist[ pos ] )
       structure_config.to( filename="POSCAR_{}.vasp".format( idx ) )


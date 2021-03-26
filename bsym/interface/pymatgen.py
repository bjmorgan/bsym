from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations, PointGroupAnalyzer
from pymatgen.util.coord import coord_list_mapping_pbc, coord_list_mapping
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure

from bsym import SpaceGroup, SymmetryOperation, ConfigurationSpace, PointGroup
from copy import copy
from functools import partial
import numpy as np

def structure_cartesian_coordinates_mapping( structure, symmop ):
    """
    Maps the coordinates of pymatgen ``Structure`` according to a ``SymmOp`` symmetry operation.

    Args:
        structure (``Structure``): The pymatgen ``Structure``.
        symmop    (``SymmOp``):    The pymatgen symmetry operation object.

    Returns
        (np.array):    The mapped Cartesian coordinates.
    """
    return structure.lattice.get_cartesian_coords( symmop.operate_multi( structure.frac_coords ) )

def molecule_cartesian_coordinates_mapping( molecule, symmop ):
    """
    Maps the coordinates of pymatgen ``Molecule`` according to a ``SymmOp`` symmetry operation.

    Args:
        molecule (``Structure``):  The pymatgen ``Molecule``.
        symmop    (``SymmOp``):    The pymatgen symmetry operation object.

    Returns
        (np.array):                The mapped Cartesian coordinates.
    """
    return symmop.operate_multi( molecule.cart_coords )

def structure_mapping_list( new_structure, mapping_structure, atol ):
    """
    Gives the index mapping between two pymatgen ``Structure`` objects.

    Args:
        new_structure (``Structure``):
        mapping_structure (``Structure``):

    Returns:
        list of indices such that mapping_structure.sites[indices] == new_structure.sites
    """
    return coord_list_mapping_pbc( new_structure.frac_coords, mapping_structure.frac_coords, atol=atol )

def molecule_mapping_list( new_molecule, mapping_molecule, atol ):
    """
    Gives the index mapping between two pymatgen ``Molecule`` objects.

    Args:
        new_structure (``Molecule``):
        mapping_structure (``Molecule``):

    Returns:
        list of indices such that mapping_molecule.sites[indices] == new_molecule.sites
    """
    return coord_list_mapping( new_molecule.cart_coords, mapping_molecule.cart_coords, atol=atol )

def unique_symmetry_operations_as_vectors_from_structure( structure, verbose=False, subset=None, atol=1e-5 ):
    """
    Uses `pymatgen`_ symmetry analysis to find the minimum complete set of symmetry operations for the space group of a structure.

    Args:
        structure (pymatgen ``Structure``): structure to be analysed.
        subset    (Optional [list]):        list of atom indices to be used for generating the symmetry operations.
        atol      (Optional [float]):       tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each symmetry operation.

    Returns:
        (list[list]): a list of lists, containing the symmetry operations as vector mappings.

    .. _pymatgen:
        http://pymatgen.org
    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc

    """
    if isinstance( structure, Structure ):
        instantiate_structure = partial( Structure, lattice=structure.lattice, coords_are_cartesian=True )
        coord_mapping = structure_cartesian_coordinates_mapping
        mapping_list = structure_mapping_list
        symmetry_analyzer = SpacegroupAnalyzer( structure )
        if verbose:
            print( "The space group for this structure is {}".format( symmetry_analyzer.get_space_group_symbol()) )
    elif isinstance( structure, Molecule ):
        instantiate_structure = Molecule
        coord_mapping = molecule_cartesian_coordinates_mapping
        mapping_list = molecule_mapping_list
        symmetry_analyzer = PointGroupAnalyzer( structure, tolerance=atol )
        if verbose:
            print( "The point group for this structure is {}".format( symmetry_analyzer.get_pointgroup()) )
    else:
        raise ValueError( 'structure argument should be a Structure or Molecule object' )
    symmetry_operations = symmetry_analyzer.get_symmetry_operations()
    mappings = []
    if subset:
        species_subset = [ spec for i,spec in enumerate( structure.species ) if i in subset ]
        cart_coords_subset = [ coord for i, coord in enumerate( structure.cart_coords ) if i in subset ]
        mapping_structure = instantiate_structure( species=species_subset, coords=cart_coords_subset )
    else:
        mapping_structure = structure
    for symmop in symmetry_operations:
        cart_coords = coord_mapping( mapping_structure, symmop )
        new_structure = instantiate_structure( species=mapping_structure.species, coords=cart_coords )
        new_mapping = [ x+1 for x in list( mapping_list( new_structure, mapping_structure, atol ) ) ]
        if new_mapping not in mappings:
            mappings.append( new_mapping )
    return mappings

def space_group_symbol_from_structure( structure ):
    """
    Returns the symbol for the space group defined by this structure. 

    Args:
        structure (pymatgen ``Structure``): The input structure.
 
    Returns:
        (str): The space group symbol.
    """
    symmetry_analyzer = SpacegroupAnalyzer( structure )
    symbol = symmetry_analyzer.get_space_group_symbol()
    return symbol

def space_group_from_structure( structure, subset=None, atol=1e-5 ):
    """
    Generates a ``SpaceGroup`` object from a `pymatgen` ``Structure``. 

    Args:
        structure (pymatgen ``Structure``): structure to be used to define the :any:`SpaceGroup`.
        subset    (Optional [list]):        list of atom indices to be used for generating the symmetry operations.
        atol      (Optional [float]):       tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each symmetry operation.

    Returns:
        a new :any:`SpaceGroup` instance

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc

    """
    mappings = unique_symmetry_operations_as_vectors_from_structure( structure, subset=subset, atol=atol )
    symmetry_operations = [ SymmetryOperation.from_vector( m ) for m in mappings ]
    return SpaceGroup( symmetry_operations=symmetry_operations )

def point_group_from_molecule( molecule, subset=None, atol=1e-5 ):
    """
    Generates a ``PointGroup`` object from a `pymatgen` ``Molecule``. 

    Args:
        molecule  (pymatgen ``Molecule``):  molecule to be used to define the :any:`PointGroup`.
        subset    (Optional [list]):        list of atom indices to be used for generating the symmetry operations.
        atol      (Optional [float]):       tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each symmetry operation.

    Returns:
        a new :any:`PointGroup` instance

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping

    """
    molecule = Molecule( molecule.species, molecule.cart_coords - molecule.center_of_mass )
    mappings = unique_symmetry_operations_as_vectors_from_structure( molecule, subset=subset, atol=atol )
    symmetry_operations = [ SymmetryOperation.from_vector( m ) for m in mappings ]
    return PointGroup( symmetry_operations=symmetry_operations )

def configuration_space_from_structure( structure, subset=None, atol=1e-5 ):
    """
    Generate a ``ConfigurationSpace`` object from a `pymatgen` ``Structure``.

    Args:
        structure (pymatgen ``Structure``): structure to be used to define the :any:`ConfigurationSpace`.
        subset    (Optional [list]):        list of atom indices to be used for generating the configuration space.
        atol      (Optional [float]):       tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each symmetry operation.
        
    Returns:
        a new :any:`ConfigurationSpace` instance.

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc

    """
    space_group = space_group_from_structure( structure, subset=subset, atol=atol )
    if subset is None:
        subset = list( range( 1, len( structure )+1 ) )
    config_space = ConfigurationSpace( objects=subset, symmetry_group=space_group )
    return config_space

def configuration_space_from_molecule( molecule, subset=None, atol=1e-5 ):
    """
    Generate a ``ConfigurationSpace`` object from a `pymatgen` ``Molecule``.

    Args:
        molecule  (pymatgen ``Molecule``):  molecule to be used to define the :any:`ConfigurationSpace`.
        subset    (Optional [list]):        list of atom indices to be used for generating the configuration space.
        atol      (Optional [float]):       tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each symmetry operation.
        
    Returns:
        a new :any:`ConfigurationSpace` instance.

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping

    """
    molecule = Molecule( molecule.species, molecule.cart_coords - molecule.center_of_mass )
    point_group = point_group_from_molecule( molecule, subset=subset, atol=atol )
    if subset is None:
        subset = list( range( 1, len( molecule )+1 ) )
    config_space = ConfigurationSpace( objects=subset, symmetry_group=point_group )
    return config_space
 
def unique_structure_substitutions( structure, to_substitute, site_distribution, verbose=False, atol=1e-5, show_progress=False ):
    """
    Generate all symmetry-unique structures formed by substituting a set of sites in a `pymatgen` structure.

    Args:
        structure (pymatgen.Structure): The parent structure.
        to_substitute (str): atom label for the sites to be substituted.
        site_distribution (dict): A dictionary that defines the number of each substituting element.
        verbose (bool): verbose output.
        atol      (Optional [float]):       tolerance factor for the ``pymatgen`` `coordinate mapping`_ under each symmetry operation. Default=1e-5.
        show_progress (opt:default=False): Show a progress bar.
                                           Setting to `True` gives a simple progress bar.
                                           Setting to `"notebook"` gives a Jupyter notebook compatible progress bar.


    Returns:
        (list[Structure]): A list of Structure objects for each unique substitution.
    
    Notes:
        The number of symmetry-equivalent configurations for each structure 
        is stored in the `number_of_equivalent_configurations` attribute. 
     
        If the parent structure was previously generated using this function
        (as part of a sequence of substitutions) the full configuration
        degeneracy of each symmetry inequivalent configuration is stored in
        the `full_configuration_degeneracy` attribute. If the parent structure
        is a standard Pymatgen Structure object, `number_of_equivalent_configurations`
        and `full_configuration_degeneracy` will be equal.

    .. _coordinate mapping:
        http://pymatgen.org/pymatgen.util.coord_utils.html#pymatgen.util.coord_utils.coord_list_mapping_pbc

    """
    site_substitution_index = list( structure.indices_from_symbol( to_substitute ) )
    if len( site_substitution_index ) != sum( site_distribution.values() ):
        raise ValueError( "Number of sites from index does not match number from site distribution" )
    if isinstance( structure, Structure ):
        config_space = configuration_space_from_structure( structure, subset=site_substitution_index, atol=atol )
    elif isinstance( structure, Molecule ):
        structure = Molecule( structure.species, structure.cart_coords - structure.center_of_mass )
        config_space = configuration_space_from_molecule( structure, subset=site_substitution_index, atol=atol )
    else:
        raise ValueError( "pymatgen Structure or Molecule object expected" )
    numeric_site_distribution, numeric_site_mapping = parse_site_distribution( site_distribution )
    unique_configurations = config_space.unique_configurations( numeric_site_distribution, verbose=verbose, show_progress=show_progress )
    new_structures = [ new_structure_from_substitution( structure, site_substitution_index, [ numeric_site_mapping[k] for k in c.tolist() ] ) for c in unique_configurations ]
    if hasattr( structure, 'number_of_equivalent_configurations' ):
        for s, c in zip( new_structures, unique_configurations ):
            s.number_of_equivalent_configurations = c.count
            s.full_configuration_degeneracy = c.count * structure.full_configuration_degeneracy
    else:
        for s, c in zip( new_structures, unique_configurations ):
            s.number_of_equivalent_configurations = c.count
            s.full_configuration_degeneracy = c.count
    return new_structures

def parse_site_distribution( site_distribution ):
    """
    Converts a site distribution using species labels into one using integer labels.

    Args:
        site_distribution (dict): e.g. `{ 'Mg': 1, 'Li': 3 }`

    Returns:
        numeric_site_distribution ( dict): e.g. `{ 1:1, 0:3 }`
        numeric_site_mapping (dict): e.g. `{ 0:'Mg', 1:'Li' }`
    """
    numeric_site_distribution = {}
    numeric_site_mapping = {}
    for i,k in enumerate( site_distribution.keys() ):
        numeric_site_distribution[i] = site_distribution[k]
        numeric_site_mapping[i] = k
    return numeric_site_distribution, numeric_site_mapping
   
def new_structure_from_substitution( parent_structure, site_substitution_index, new_species_list ):
    """
    Generate a new pymatgen ``Structure`` from site substitution parameters.

    Args:
        parent_structure (Structure):        The parent pymatgen ``Struture`` object.
        site_substitution_index (list[int]): The list of site indices to be substituted.
        new_species_list (list[str]):        A list of the replacement atomic species.
    
    Returns:
        (``Structure``): The new pymatgen ``Structure``.

    Notes:
        pymatgen ``Structure`` and ``Molecule`` classes both subclass ``SiteCollection``.
        This function will also accept a parent ``Molecule`` object, and return a new
        ``Molecule``.
    """
    if len( site_substitution_index ) != len( new_species_list ):
        raise ValueError
    if any( i >= len( parent_structure ) for i in site_substitution_index ):
        raise ValueError
    s = parent_structure.copy()
    for i, spec in zip( site_substitution_index, new_species_list ):
        s[i] = spec
    return s


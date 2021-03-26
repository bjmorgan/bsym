import unittest
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.operations import SymmOp
from bsym.interface.pymatgen import ( unique_symmetry_operations_as_vectors_from_structure, 
                                      space_group_from_structure, 
                                      parse_site_distribution, 
                                      unique_structure_substitutions, 
                                      new_structure_from_substitution, 
                                      configuration_space_from_structure, 
                                      space_group_symbol_from_structure, 
                                      configuration_space_from_molecule, 
                                      structure_cartesian_coordinates_mapping,
                                      molecule_cartesian_coordinates_mapping )

from itertools import permutations
from bsym import SymmetryOperation, Configuration, SpaceGroup, PointGroup, ConfigurationSpace

class TestPymatgenInterface( unittest.TestCase ):

    def setUp( self ):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.0, 0.5, 0.5 ],
                             [ 0.5, 0.0, 0.5 ] ] )
        atom_list = [ 'Li' ] * len( coords )
        lattice = Lattice.from_parameters( a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        self.structure = Structure( lattice, atom_list, coords )
        # construct a pymatgen Molecule instance
        # square molecule (D4h)
        m_coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                               [ 1.0, 0.0, 0.0 ],
                               [ 0.0, 1.0, 0.0 ],
                               [ 1.0, 1.0, 0.0 ] ] )
        molecule = Molecule( atom_list, m_coords )
        molecule = Molecule( molecule.species, molecule.cart_coords - molecule.center_of_mass )
        self.molecule = molecule 

    def test_new_structure_from_substitution( self ):
        substitution_index = [ 2,3 ]
        new_species_list = [ 'Mg', 'Fe' ] 
        s_new = new_structure_from_substitution( self.structure, substitution_index, new_species_list ) 
        self.assertEqual( s_new[2].species_string, 'Mg' )
        self.assertEqual( s_new[3].species_string, 'Fe' )

    def test_new_structure_from_substitution_raises_ValueError_with_oversize_index( self ):
        substitution_index = [ 0, 1, 2, 3, 4 ]
        new_species_list = [ 'Mg', 'Fe' ]
        with self.assertRaises( ValueError ):
            new_structure_from_substitution( self.structure, substitution_index, new_species_list )

    def test_new_structure_from_substitution_raises_ValueError_with_invalid_index( self ):
        substitution_index = [ 2, 4 ]
        new_species_list = [ 'Mg', 'Fe' ]
        with self.assertRaises( ValueError ):
            new_structure_from_substitution( self.structure, substitution_index, new_species_list )

    def test_parse_site_distribution( self ):
        site_distribution = { 'Mg': 1, 'Li': 3 }
        n, d = parse_site_distribution( site_distribution )
        for k, v in n.items():
            self.assertEqual( site_distribution[ d[ k ] ], v )

    def test_structure_cartesian_coordinates_mapping( self ):
        mock_symmop = Mock( spec=SymmOp )
        new_coords = np.array( [ [ 0.5, 0.5, 0.5 ] ] )
        mock_symmop.operate_multi = Mock( return_value=new_coords )
        self.structure.lattice.get_cartesian_coords = Mock( return_value=np.array( [ [ 2.0, 2.0, 2.0 ] ] ) )
        mapped_coords = structure_cartesian_coordinates_mapping( self.structure, mock_symmop )
        np.testing.assert_array_equal( mapped_coords, np.array( [ [ 2.0, 2.0, 2.0 ] ] ) )
        np.testing.assert_array_equal( mock_symmop.operate_multi.call_args[0][0], self.structure.frac_coords )

    def test_molecule_cartesian_coordinates_mapping( self ):
        mock_symmop = Mock( spec=SymmOp )
        new_coords = np.array( [ [ 0.5, 0.5, 0,5 ] ] )
        mock_symmop.operate_multi = Mock( return_value=new_coords )
        mapped_coords = molecule_cartesian_coordinates_mapping( self.molecule, mock_symmop ) 
        np.testing.assert_array_equal( mapped_coords, new_coords )
        np.testing.assert_array_equal( mock_symmop.operate_multi.call_args[0][0], self.molecule.cart_coords )

if __name__ == '__main__':
    unittest.main()

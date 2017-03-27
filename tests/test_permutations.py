import unittest
from bsym import permutations

class PermutationsTestCase( unittest.TestCase ):
    """Tests for permutations functions"""

    def test_flatten_list( self ):
        l = [ [ 1, 2 ], [ 3, 4, 5 ] ]
        self.assertEqual( permutations.flatten_list( l ), [ 1, 2, 3, 4, 5 ] )

    def test_unique_permuations( self ):
        seq = [ 1, 1, 0, 0 ]
        all_permutations = [ [ 1, 1, 0, 0 ],
                             [ 1, 0, 1, 0 ],
                             [ 1, 0, 0, 1 ],
                             [ 0, 1, 1, 0 ],
                             [ 0, 1, 0, 1 ],
                             [ 0, 0, 1, 1 ] ]
        unique_permutations = list( permutations.unique_permutations( seq ) )
        for p in all_permutations:
            self.assertEqual( p in all_permutations, True )
        self.assertEqual( len( all_permutations ), len( unique_permutations ) )
 
if __name__ == '__main__':
    unittest.main()

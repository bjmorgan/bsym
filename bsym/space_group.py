from bsym import SymmetryGroup

class SpaceGroup( SymmetryGroup ):

    def __repr__( self ):
        to_return = 'SpaceGroup\n'
        for so in self.symmetry_operations:
            to_return += f"{so.label}\t{so.as_vector()}\n"
        return to_return

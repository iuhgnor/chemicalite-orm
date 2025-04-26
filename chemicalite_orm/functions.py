from sqlalchemy.dialects import sqlite
from sqlalchemy.sql.functions import GenericFunction
from sqlalchemy.types import Float, Integer

from .types import Bfp, Mol


class mol_from_smarts(GenericFunction):
    inherit_cache = True

    name = "mol_from_smarts"
    type = Mol()


class mol_is_substruct(GenericFunction):
    inherit_cache = True

    name = "mol_is_substruct"
    type = Integer()


class mol_compare(GenericFunction):
    name = "mol_compare"
    type = Integer()


class bfp_tanimoto(GenericFunction):
    name = "bfp_tanimoto"
    type = Float()


class is_substructure(GenericFunction):
    name = "is_substructure"
    type = Integer()


class tanimoto(GenericFunction):
    name = "tanimoto"
    type = Float()


class mol_pattern_bfp(GenericFunction):
    type = Bfp()
    name = "mol_pattern_bfp"


class rdtree_subset(GenericFunction):
    type = Bfp()
    name = "rdtree_subset"


class mol_from_binary_mol(GenericFunction):
    name = "mol_from_binary_mol"
    type = Mol()


class mol_to_binary_mol(GenericFunction):
    name = "mol_to_binary_mol"
    type = sqlite.BLOB()


class mol_from_smiles(GenericFunction):
    inherit_cache = True

    name = "mol_from_smiles"
    type = Mol()

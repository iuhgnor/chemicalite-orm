from sqlalchemy.sql.functions import GenericFunction
from sqlalchemy.types import Float, Integer

from .types import Mol


class mol_from_smarts(GenericFunction):
    inherit_cache = True

    name = "mol_from_smarts"
    type = Mol()


class mol_is_substruct(GenericFunction):
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

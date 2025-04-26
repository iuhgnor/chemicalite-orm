from rdkit import Chem
from sqlalchemy.dialects import sqlite
from sqlalchemy.sql.functions import GenericFunction
from sqlalchemy.types import UserDefinedType


class Mol(UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "mol"

    def bind_processor(self, dialect):
        def process(value):
            if isinstance(value, Chem.Mol):
                return Chem.MolToSmiles(value)
            elif isinstance(value, str):
                return value
            return None

        return process

    def bind_expression(self, bindvalue):
        return mol_from_smiles(bindvalue)

    def column_expression(self, col):
        return mol_to_binary_mol(col)

    def result_processor(self, dialect, coltype):
        def process(value):
            if value is None:
                return None
            return Chem.Mol(value)

        return process


class Bfp(UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "bfp"


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

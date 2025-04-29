from rdkit import Chem
from sqlalchemy import func
from sqlalchemy.types import UserDefinedType


class Mol(UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "mol"

    def bind_processor(self, dialect):
        def process(value):
            if isinstance(value, Chem.Mol):
                return value.ToBinary()
            elif isinstance(value, str):
                return Chem.MolFromSmiles(value).ToBinary()
            return None

        return process

    def bind_expression(self, bindvalue):
        return func.mol_from_binary_mol(bindvalue, type_=self)

    def column_expression(self, col):
        # ref https://groups.google.com/g/sqlalchemy/c/2nGwTuwkNxw?pli=1
        return func.mol_to_binary_mol(col, type_=self)

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

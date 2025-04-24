from typing import Optional

from pydantic import ConfigDict
from rdkit import Chem, DataStructs
from rdkit.DataStructs import ExplicitBitVect
from sqlalchemy.sql.functions import GenericFunction
from sqlalchemy.types import Float, Integer, UserDefinedType
from sqlmodel import Field, SQLModel


class Mol(UserDefinedType):
    cache_ok = True

    def hassubstruct(self, other):
        return mol_is_substruct(other, self.expr) == 1

    def issubstruct(self, other):
        return mol_is_substruct(self.expr, other) == 1

    def compare(self, other):
        return mol_compare(self.expr, other) == 1

    def get_col_spec(self, **kw):
        return "mol"

    def bind_processor(self, dialect):
        def process(value):
            if isinstance(value, Chem.Mol):
                value = memoryview(value.ToBinary())
            elif isinstance(value, str):
                value = memoryview(Chem.MolFromSmiles(value).ToBinary())
            return value

        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            if value is None:
                return value
            return Chem.Mol(bytes(value))

        return process


class Bfp(UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "bfp"

    def bind_processor(self, dialect):
        def process(value):
            if isinstance(value, ExplicitBitVect):
                return DataStructs.BitVectToBinaryText(value)
            return value

        return process

    def result_processor(self, dialect, coltype):
        def process(value):
            if value is None:
                return None
            return DataStructs.CreateFromBinaryText(value)

        return process


class mol_from_smiles(GenericFunction):
    name = "mol_from_smiles"
    type = Mol()


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


class mol_morgan_bfp(GenericFunction):
    inherit_cache = True

    name = "mol_morgan_bfp"
    type = Bfp()


class bfp_tanimoto(GenericFunction):
    name = "bfp_tanimoto"
    type = Float()


class is_substructure(GenericFunction):
    name = "is_substructure"
    type = Integer()


class tanimoto(GenericFunction):
    name = "tanimoto"
    type = Float()


class Compound(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: str
    molecule: Mol = Field(sa_type=Mol)
    mfp2: Bfp = Field(sa_type=Bfp)

    model_config = ConfigDict(arbitrary_types_allowed=True)  # type: ignore

    def __repr__(self):
        return f"({self.name}) < {self.smiles} >"

    def is_substructure_of(self, other_mol: Mol) -> bool:
        return bool(mol_is_substruct(other_mol, self.molecule))

    def tanimoto_similarity(self, other_mol: Mol) -> float:
        return float(bfp_tanimoto(self.mfp2, mol_morgan_bfp(other_mol)))

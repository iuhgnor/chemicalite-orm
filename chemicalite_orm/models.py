from typing import Optional

from pydantic import ConfigDict
from sqlmodel import Field, SQLModel

from .types import Mol


class Compound(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: str
    molecule: Mol = Field(sa_type=Mol)

    model_config = ConfigDict(arbitrary_types_allowed=True)  # type: ignore

    def __repr__(self):
        return f"({self.name}) < {self.smiles} >"

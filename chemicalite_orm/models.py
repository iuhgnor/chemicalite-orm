from typing import Optional, Type

from pydantic import ConfigDict
from sqlalchemy import event, text
from sqlalchemy.orm import Session as ORMSession
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


# === Auto Index Utilities ===
def create_virtual_index_for_model(conn, model: Type[SQLModel]):
    table_name = str(model.__tablename__)
    index_name = f"str_idx_{table_name}"
    conn.execute(
        text(f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {index_name}
        USING rdtree(id, fp bits(2048))
    """)
    )
    # conn.execute(
    #     text(f"""
    #     INSERT INTO {index_name}(id, fp)
    #     SELECT id, mol_pattern_bfp(molecule, {bits})
    #     FROM {table_name} WHERE molecule IS NOT NULL
    # """)
    # )


def register_auto_index(model: Type[SQLModel]):
    @event.listens_for(model.__table__, "after_create")  # type: ignore
    def auto_index_creation(target, connection, **kw):
        create_virtual_index_for_model(connection, model)

    @event.listens_for(ORMSession, "after_flush")
    def auto_update_index(session, flush_context):
        index_name = f"str_idx_{model.__tablename__}"
        for instance in session.new:
            if isinstance(instance, model):
                stmt = text(f"""
                    INSERT INTO {index_name}(id, fp)
                    VALUES (:id, mol_pattern_bfp(mol_from_smiles(:mol), 2048))
                """)
                session.execute(stmt, {"id": instance.id, "mol": instance.molecule})
        for instance in session.deleted:
            if isinstance(instance, model):
                stmt = text(f"DELETE FROM {index_name} WHERE id = :id")
                session.execute(stmt, {"id": instance.id})


# Register auto index for Compound
register_auto_index(Compound)

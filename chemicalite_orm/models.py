from typing import Optional, Type

from pydantic import ConfigDict
from sqlalchemy import event, text
from sqlalchemy.orm import Session as ORMSession
from sqlmodel import Field, Session, SQLModel

from .types import Mol


class Compound(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: str
    molecule: Mol = Field(sa_type=Mol)

    model_config = ConfigDict(arbitrary_types_allowed=True)  # type: ignore

    def __repr__(self):
        return f"({self.name}) < {self.smiles} >"

    @classmethod
    def search_by_substructure(cls, session: Session, submol_smiles: str):
        """Substructure search using RDKit fingerprint."""
        table_name = str(cls.__tablename__)
        idx_tabel_name = f"str_idx_{table_name}"
        sql = text(f"""
            SELECT {table_name}.* FROM {table_name}
            JOIN {idx_tabel_name} ON compound.id = {idx_tabel_name}.id
            WHERE {idx_tabel_name}.id MATCH rdtree_subset(mol_pattern_bfp(mol_from_smiles('{submol_smiles}'), 2048));
            """)
        return session.exec(sql).all()  # type: ignore

    @classmethod
    def search_by_similarity(
        cls, session: Session, query_smiles: str, threshold: float = 0.7
    ):
        """Similarity search using Tanimoto coefficient."""
        table_name = str(cls.__tablename__)
        sim_table_name = f"sim_idx_{table_name}"
        sql = text(f"""
            SELECT {table_name}.*, bfp_tanimoto({sim_table_name}.fp, mol_morgan_bfp(mol_from_smiles('{query_smiles}'), 2, 1024)) as sim
            FROM {table_name}
            JOIN {sim_table_name} ON {table_name}.id = {sim_table_name}.id
            WHERE sim >= {threshold}
            ORDER BY sim DESC
        """)

        return session.exec(sql).all()  # type: ignore


# === Auto Index Utilities ===
def create_virtual_index_for_model(conn, model: Type[SQLModel]):
    table_name = str(model.__tablename__)
    index_name = f"str_idx_{table_name}"
    sql = text(f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {index_name}
        USING rdtree(id, fp bits(2048))
    """)
    conn.execute(sql)  # type: ignore

    sim_index_name = f"sim_idx_{table_name}"
    sql = text(f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {sim_index_name}
        USING rdtree(id, fp bits(1024))
    """)
    conn.execute(sql)  # type: ignore


def register_auto_index(model: Type[SQLModel]):
    @event.listens_for(model.__table__, "after_create")  # type: ignore
    def auto_index_creation(target, connection, **kw):
        create_virtual_index_for_model(connection, model)

    @event.listens_for(ORMSession, "after_flush")
    def auto_update_index(session, flush_context):
        table_name = model.__tablename__
        str_index_name = f"str_idx_{table_name}"
        sim_index_name = f"sim_idx_{table_name}"
        for instance in session.new:
            if isinstance(instance, model):
                stmt1 = text(f"""
                    INSERT INTO {str_index_name}(id, fp)
                    VALUES (:id, mol_pattern_bfp(mol_from_smiles(:mol), 2048))
                """)
                stmt2 = text(f"""
                    INSERT INTO {sim_index_name}(id, fp)
                    VALUES (:id, mol_morgan_bfp(mol_from_smiles(:mol), 2, 1024))
                """)
                session.execute(stmt1, {"id": instance.id, "mol": instance.smiles})
                session.execute(stmt2, {"id": instance.id, "mol": instance.smiles})
        for instance in session.deleted:
            if isinstance(instance, model):
                session.execute(
                    text(f"DELETE FROM {str_index_name} WHERE id = :id"),
                    {"id": instance.id},
                )
                session.execute(
                    text(f"DELETE FROM {sim_index_name} WHERE id = :id"),
                    {"id": instance.id},
                )


# Register auto index for Compound
register_auto_index(Compound)

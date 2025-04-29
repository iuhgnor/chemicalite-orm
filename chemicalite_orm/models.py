from typing import Type

from sqlalchemy import Integer, String, event, text
from sqlalchemy.orm import Mapped, Session, declarative_base, mapped_column

from .types import Mol

Base = declarative_base()


class Compound(Base):  # type: ignore
    __tablename__ = "compound"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String, nullable=False)
    smiles: Mapped[str] = mapped_column(String, nullable=False)
    molecule: Mapped[Mol] = mapped_column(Mol, nullable=False)

    def __repr__(self):
        return f"({self.name}) < {self.smiles} >"

    @classmethod
    def search_by_substructure(cls, session: Session, submol_smiles: str):
        table_name = cls.__tablename__
        idx_table_name = f"str_idx_{table_name}"
        sql = text(f"""
            SELECT {table_name}.* FROM {table_name}
            JOIN {idx_table_name} ON {table_name}.id = {idx_table_name}.id
            WHERE {idx_table_name}.id MATCH rdtree_subset(mol_pattern_bfp(mol_from_smiles(:querymol), 2048));
        """)
        return session.execute(sql, {"querymol": submol_smiles}).fetchall()

    @classmethod
    def search_by_similarity(
        cls, session: Session, query_smiles: str, threshold: float = 0.7
    ):
        table_name = cls.__tablename__
        sim_table_name = f"sim_idx_{table_name}"
        sql = text(f"""
            SELECT {table_name}.*, bfp_tanimoto({sim_table_name}.fp, mol_morgan_bfp(mol_from_smiles(:querymol), 2, 1024)) as sim
            FROM {table_name}
            JOIN {sim_table_name} ON {table_name}.id = {sim_table_name}.id
            WHERE sim >= :threshold
            ORDER BY sim DESC
        """)
        return session.execute(
            sql, {"querymol": query_smiles, "threshold": threshold}
        ).fetchall()


# === Auto Index Utilities ===
def create_virtual_index_for_model(conn, model: Type[Base]):  # type: ignore
    table_name = model.__tablename__  # type: ignore
    str_index_name = f"str_idx_{table_name}"
    sim_index_name = f"sim_idx_{table_name}"

    conn.execute(
        text(f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {str_index_name}
        USING rdtree(id, fp bits(2048))
    """)
    )

    conn.execute(
        text(f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {sim_index_name}
        USING rdtree(id, fp bits(1024))
    """)
    )


def register_auto_index(model: Type[Base]):  # type: ignore
    @event.listens_for(model.__table__, "after_create")  # type: ignore
    def auto_index_creation(target, connection, **kw):
        create_virtual_index_for_model(connection, model)

    @event.listens_for(Session, "after_flush")
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

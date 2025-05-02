from sqlalchemy import Integer, String, text
from sqlalchemy.orm import Mapped, Session, declarative_base, mapped_column

from .search import enable_chem_search
from .types import Mol

Base = declarative_base()


@enable_chem_search()
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
            WHERE {idx_table_name}.id MATCH rdtree_subset(mol_pattern_bfp(mol_from_smiles('{submol_smiles}'), 2048));
        """)
        return session.execute(sql).fetchall()

    @classmethod
    def search_by_similarity(
        cls, session: Session, query_smiles: str, threshold: float = 0.7
    ):
        table_name = cls.__tablename__
        sim_table_name = f"sim_idx_{table_name}"
        sql = text(f"""
            SELECT {table_name}.*, bfp_tanimoto({sim_table_name}.fp, mol_morgan_bfp(mol_from_smiles('{query_smiles}'), 2, 2048)) as sim
            FROM {table_name}
            JOIN {sim_table_name} ON {table_name}.id = {sim_table_name}.id
            WHERE sim >= {threshold}
            ORDER BY sim DESC
        """)
        return session.execute(sql).fetchall()

from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy.orm import DeclarativeMeta
from sqlalchemy.sql import text


def enable_chem_search(fp_bits: int = 2048, fp_radius: int = 2):
    def decorator(cls: DeclarativeMeta):
        table_name = cls.__tablename__  # type: ignore
        str_idx_table = f"str_idx_{table_name}"
        sim_idx_table = f"sim_idx_{table_name}"

        @hybrid_method
        def substructure_search(cls, submol_smiles: str):
            sql = f"""
                {table_name}.id IN (
                    SELECT id FROM {str_idx_table}
                    WHERE id MATCH rdtree_subset(
                        mol_pattern_bfp(mol_from_smiles(:querymol), {fp_bits})
                    )
                )
            """
            return text(sql).bindparams(querymol=submol_smiles)

        @hybrid_method
        def similarity_search(cls, query_smiles: str, threshold: float = 0.7):
            sql = f"""
                {table_name}.id IN (
                    SELECT id FROM {sim_idx_table}
                    WHERE bfp_tanimoto(
                        fp,
                        mol_morgan_bfp(mol_from_smiles(:querymol), {fp_radius}, {fp_bits})
                    ) >= :threshold
                )
            """
            return text(sql).bindparams(querymol=query_smiles, threshold=threshold)

        cls.substructure_search = substructure_search
        cls.similarity_search = similarity_search

        return cls

    return decorator

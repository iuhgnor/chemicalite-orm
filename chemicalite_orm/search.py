from sqlalchemy import event
from sqlalchemy.ext.hybrid import hybrid_method
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.sql import text


def enable_chem_search(
    mol_col: str = "molecule",
    id_col: str = "id",
    fp_bits: int = 2048,
    fp_radius: int = 2,
):
    def decorator(cls: DeclarativeBase):
        table_name = cls.__tablename__  # type: ignore
        str_idx_table = f"str_idx_{table_name}"
        sim_idx_table = f"sim_idx_{table_name}"

        @hybrid_method
        def has_substructure(cls, submol_smiles: str):
            sql = f"""
                {table_name}.{id_col} IN (
                    SELECT id FROM {str_idx_table}
                    WHERE id MATCH rdtree_subset(
                        mol_pattern_bfp(mol_from_smiles(:querymol), {fp_bits})
                    )
                )
            """
            return text(sql).bindparams(querymol=submol_smiles)

        @hybrid_method
        def similarity_to(cls, query_smiles: str, threshold: float = 0.7):
            sql = f"""
                {table_name}.{id_col} IN (
                    SELECT id FROM {sim_idx_table}
                    WHERE bfp_tanimoto(
                        fp,
                        mol_morgan_bfp(mol_from_smiles(:querymol), {fp_radius}, {fp_bits})
                    ) >= :threshold
                )
            """
            return text(sql).bindparams(querymol=query_smiles, threshold=threshold)

        cls.has_substructure = has_substructure  # type: ignore
        cls.similarity_to = similarity_to  # type: ignore

        @event.listens_for(cls.__table__, "after_create")
        def create_index_tables(target, connection, **kw):
            connection.execute(
                text(f"""
                CREATE VIRTUAL TABLE IF NOT EXISTS {str_idx_table}
                USING rdtree(id, fp bits({fp_bits}))
            """)
            )

            connection.execute(
                text(f"""
                CREATE VIRTUAL TABLE IF NOT EXISTS {sim_idx_table}
                USING rdtree(id, fp bits({fp_bits}))
            """)
            )

            connection.execute(
                text(f"""
                CREATE TRIGGER IF NOT EXISTS trg_insert_{table_name}
                AFTER INSERT ON {table_name}
                BEGIN
                    INSERT INTO {str_idx_table}(id, fp)
                    VALUES (NEW.{id_col}, mol_pattern_bfp(NEW.{mol_col}, {fp_bits}));

                    INSERT INTO {sim_idx_table}(id, fp)
                    VALUES (NEW.{id_col}, mol_morgan_bfp(NEW.{mol_col}, {fp_radius}, {fp_bits}));
                END;
            """)
            )

            connection.execute(
                text(f"""
                CREATE TRIGGER IF NOT EXISTS trg_delete_{table_name}
                AFTER DELETE ON {table_name}
                BEGIN
                    DELETE FROM {str_idx_table} WHERE id = OLD.{id_col};
                    DELETE FROM {sim_idx_table} WHERE id = OLD.{id_col};
                END;
            """)
            )

            connection.execute(
                text(f"""
                CREATE TRIGGER IF NOT EXISTS trg_update_{table_name}
                AFTER UPDATE OF {mol_col} ON {table_name}
                BEGIN
                    DELETE FROM {str_idx_table} WHERE id = OLD.{id_col};
                    DELETE FROM {sim_idx_table} WHERE id = OLD.{id_col};

                    INSERT INTO {str_idx_table}(id, fp)
                    VALUES (NEW.{id_col}, mol_pattern_bfp(NEW.{mol_col}, {fp_bits}));

                    INSERT INTO {sim_idx_table}(id, fp)
                    VALUES (NEW.{id_col}, mol_morgan_bfp(NEW.{mol_col}, {fp_radius}, {fp_bits}));
                END;
            """)
            )

        return cls

    return decorator

from sqlalchemy import text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import DeclarativeBase


def setup_mol_index_triggers(
    engine: Engine,
    model: DeclarativeBase,
    id_column: str = "id",
    mol_column: str = "molecule",
    fp_bits: int = 2048,
    fp_radius: int = 2,
):
    table_name = model.__tablename__
    str_idx = f"str_idx_{table_name}"
    sim_idx = f"sim_idx_{table_name}"

    statements = [
        f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {str_idx}
        USING rdtree(id, fp bits({fp_bits}));
        """,
        f"""
        CREATE VIRTUAL TABLE IF NOT EXISTS {sim_idx}
        USING rdtree(id, fp bits({fp_bits}));
        """,
        f"""
        CREATE TRIGGER IF NOT EXISTS insert_{table_name}_index
        AFTER INSERT ON {table_name}
        BEGIN
          INSERT INTO {str_idx}(id, fp)
          VALUES (
            NEW.{id_column},
            mol_pattern_bfp(NEW.{mol_column}, {fp_bits})
          );

          INSERT INTO {sim_idx}(id, fp)
          VALUES (
            NEW.{id_column},
            mol_morgan_bfp(NEW.{mol_column}, {fp_radius}, {fp_bits})
          );
        END;
        """,
        f"""
        CREATE TRIGGER IF NOT EXISTS delete_{table_name}_index
        AFTER DELETE ON {table_name}
        BEGIN
          DELETE FROM {str_idx} WHERE id = OLD.{id_column};
          DELETE FROM {sim_idx} WHERE id = OLD.{id_column};
        END;
        """,
        f"""
        CREATE TRIGGER IF NOT EXISTS update_{table_name}_index
        AFTER UPDATE OF {mol_column} ON {table_name}
        BEGIN
          DELETE FROM {str_idx} WHERE id = OLD.{id_column};
          DELETE FROM {sim_idx} WHERE id = OLD.{id_column};

          INSERT INTO {str_idx}(id, fp)
          VALUES (
            NEW.{id_column},
            mol_pattern_bfp(NEW.{mol_column}, {fp_bits})
          );

          INSERT INTO {sim_idx}(id, fp)
          VALUES (
            NEW.{id_column},
            mol_morgan_bfp(NEW.{mol_column}, {fp_radius}, {fp_bits})
          );
        END;
        """,
    ]

    with engine.begin() as conn:
        for statement in statements:
            conn.execute(text(statement))

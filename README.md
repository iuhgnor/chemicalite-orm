# ChemicaLite ORM

基于[raiz](https://github.com/rvianello/razi)和[sqlmodel-rdkit-example](https://github.com/savvan0h/sqlmodel-rdkit-example)代码, 开发SQLite数据库RDKit cartridge (ChemicaLite) ORM.

## 安装

```python
conda create -n test python=3.9
conda activate test

conda install rdkit chemicalite -c conda-forge
pip install git+https://github.com/iuhgnor/chemicalite-orm.git
```

## 使用

1. 创建数据库

    ```python
    from sqlalchemy import Integer, String
    from sqlalchemy.orm import Mapped, declarative_base, mapped_column

    from chemicalite_orm.search import enable_chem_search
    from chemicalite_orm.types import Mol

    DATABASE_URL = "sqlite:///chemicalite.db"

    engine = create_engine(DATABASE_URL, echo=True)


    @event.listens_for(engine, "connect")
    def load_chemicalite(dbapi_conn, connection_record):
        dbapi_conn.enable_load_extension(True)
        dbapi_conn.load_extension("chemicalite")


    Base = declarative_base()


    @enable_chem_search(mol_col="molecule", id_col="id", fp_bits=2048, fp_radius=2)
    class Compound(Base):  # type: ignore
        __tablename__ = "compound"

        id: Mapped[int] = mapped_column(Integer, primary_key=True)
        name: Mapped[str] = mapped_column(String, nullable=False)
        smiles: Mapped[str] = mapped_column(String, nullable=False)
        molecule: Mapped[Mol] = mapped_column(Mol, nullable=False)


    Base.metadata.create_all(engine)       
    SessionLocal = sessionmaker(bind=engine)

    ```

1. 插入数据

    ```python
    with SessionLocal() as session:
        for i, smiles in enumerate(mols):
            compound = Compound(
                name=f"mol_{i}",
                smiles=smiles,
                molecule=Chem.MolFromSmiles(smiles),
            )
            session.add(compound)
        session.commit()

    ```

2. 检索数据

   ```python
    with SessionLocal() as session:
        stmt = select(Compound).where(Compound.name == "mol_1")
        result = session.execute(stmt).scalars().first()
        mol = result.molecule

    mol
   ```

   ![output](./examples/images/3.png)

3. 子结构查询

    ```python
    with SessionLocal() as session:
        stmt = select(Compound).where(Compound.has_substructure("c1ccccc1"))
        results = session.execute(stmt).scalars().all()

    mols = [r.molecule for r in results[:10]]
    Draw.MolsToGridImage(mols, molsPerRow=5)

    ```

    ![output](./examples/images/4.png)

5. 相似性查询

    ```python
    with SessionLocal() as session:
        stmt = select(Compound).where(Compound.similarity_to("c1ccnnc1", 0.1))
        results = session.execute(stmt).scalars().all()

    mols = [r.molecule for r in results[:10]]
    Draw.MolsToGridImage(mols, molsPerRow=5)


    ```

    ![output](./examples/images/5.png)

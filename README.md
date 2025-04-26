# ChemicaLite ORM

基于[raiz](https://github.com/rvianello/razi)和[sqlmodel-rdkit-example](https://github.com/savvan0h/sqlmodel-rdkit-example)代码, 开发SQLite数据库RDKit cartridge (ChemicaLite) ORM.

## 安装

```python
conda create -n test python=3.9
conda activate test

conda install rdkit chemicalite -c conda-forge
pip install sqlmodel
pip install git+https://github.com/savvan0h/chemicaLite-orm.git
```

## 使用

1. 创建数据库

    ```python
    from sqlmodel import create_engine
    from sqlalchemy import event

    from chemicalite_orm.models import Compound


    DATABASE_URL = "sqlite:///chemicalite.db"
    engine = create_engine(DATABASE_URL, echo=True)

    # 加载ChemicaLite
    @event.listens_for(engine, "connect")
    def load_chemicalite(dbapi_conn, connection_record):
        dbapi_conn.enable_load_extension(True)
        dbapi_conn.load_extension("chemicalite")


    Compound.__table__.create(engine)

    ```

2. 插入数据

    ```python
    from sqlmodel import Session

    from chemicalite_orm.molecules import SMILES_SAMPLE


    with Session(engine) as session:
        for i, smiles in enumerate(SMILES_SAMPLE):
            compound = Compound(
                name=f"mol_{i}",
                smiles=smiles,
                molecule=smiles,
            )
            session.add(compound)
        session.commit()

    ```

3. 子结构查询

    ```python
    from rdkit import Chem
    from rdkit.Chem import Draw


    with Session(engine) as session:
        hits = Compound.search_by_substructure(session, "c1ccnnc1")
        print(len(hits))
        

    hit_smis = [hit.smiles for hit in hits]
    hit_mols = [Chem.MolFromSmiles(smi) for smi in hit_smis[:10]]
    Draw.MolsToGridImage(hit_mols, molsPerRow=5)

    ```

    ![output](./examples/images/1.png)

4. 相似性查询

    ```python
    from rdkit import Chem
    from rdkit.Chem import Draw


    with Session(engine) as session:
        hits = Compound.search_by_similarity(session, "c1ccnnc1", 0.01)
        print(len(hits))


    hit_smis = [hit.smiles for hit in hits[:10]]
    hit_mols = [Chem.MolFromSmiles(smi) for smi in hit_smis]
    similarity = [hit.sim for hit in hits[:10]]
    Draw.MolsToGridImage(
        hit_mols, molsPerRow=5, legends=[f"{s:.2f}" for s in similarity]
    )


    ```

    ![output](./examples/images/2.png)

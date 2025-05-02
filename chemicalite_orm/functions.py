from sqlalchemy.dialects.sqlite import BLOB
from sqlalchemy.sql import functions
from sqlalchemy.types import TEXT, Float, Integer

from .types import Bfp, Mol

_FUNCTIONS = [
    ("mol_from_smiles", Mol, None),
    ("mol_from_smarts", Mol, None),
    ("mol_from_molblock", Mol, None),
    ("mol_from_binary_mol", Mol, None),
    ("mol_to_smiles", TEXT, None),
    ("mol_to_smarts", TEXT, None),
    ("mol_to_molblock", TEXT, None),
    ("mol_to_binary_mol", BLOB, None),
    ("mol_is_substruct", Integer, None),
    ("mol_is_superstruct", Integer, None),
    ("mol_cmp", Integer, None),
    ("mol_hba", Integer, None),
    ("mol_hbd", Integer, None),
    ("mol_num_atms", Integer, None),
    ("mol_num_hvyatms", Integer, None),
    ("mol_num_rotatable_bnds", Integer, None),
    ("mol_num_hetatms", Integer, None),
    ("mol_num_rings", Integer, None),
    ("mol_num_aromatic_rings", Integer, None),
    ("mol_num_aliphatic_rings", Integer, None),
    ("mol_num_saturated_rings", Integer, None),
    ("mol_mw", Float, None),
    ("mol_tpsa", Float, None),
    ("mol_fraction_csp3", Float, None),
    ("mol_chi0v", Float, None),
    ("mol_chi0n", Float, None),
    ("mol_kappa1", Float, None),
    ("mol_logp", Float, None),
    ("mol_formula", TEXT, None),
    ("mol_hash_anonymousgraph", TEXT, None),
    ("mol_hash_elementgraph", TEXT, None),
    ("mol_hash_canonicalsmiles", TEXT, None),
    ("mol_hash_murckoscaffold", TEXT, None),
    ("mol_hash_extendedmurcko", TEXT, None),
    ("mol_hash_molformula", TEXT, None),
    ("mol_hash_atombondcounts", TEXT, None),
    ("mol_hash_degreevector", TEXT, None),
    ("mol_hash_mesomer", TEXT, None),
    ("mol_hash_hetatomtautomer", TEXT, None),
    ("mol_hash_hetatomprotomer", TEXT, None),
    ("mol_hash_redoxpair", TEXT, None),
    ("mol_hash_regioisomer", TEXT, None),
    ("mol_hash_netcharge", TEXT, None),
    ("mol_hash_smallworldindexbr", TEXT, None),
    ("mol_hash_smallworldindexbrl", TEXT, None),
    ("mol_hash_arthorsubstructureorder", TEXT, None),
    ("mol_prop_list", list[TEXT], None),
    ("mol_has_prop", Integer, None),
    ("mol_set_prop", Mol, None),
    ("mol_get_TEXT_prop", TEXT, None),
    ("mol_get_Integer_prop", Integer, None),
    ("mol_get_float_prop", Float, None),
    ("mol_delete_substructs", Mol, None),
    ("mol_replace_substructs", list[Mol], None),
    ("mol_replace_sidechains", Mol, None),
    ("mol_replace_core", Mol, None),
    ("mol_murcko_decompose", Mol, None),
    ("mol_find_mcs", Mol, None),
    ("mol_cleanup", Mol, None),
    ("mol_normalize", Mol, None),
    ("mol_reionize", Mol, None),
    ("mol_remove_fragments", Mol, None),
    ("mol_canonical_tautomer", Mol, None),
    ("mol_tautomer_parent", Mol, None),
    ("mol_fragment_parent", Mol, None),
    ("mol_stereo_parent", Mol, None),
    ("mol_isotope_parent", Mol, None),
    ("mol_charge_parent", Mol, None),
    ("mol_super_parent", Mol, None),
    ("mol_layered_Bfp", Bfp, None),
    ("mol_rdkit_Bfp", Bfp, None),
    ("mol_atom_pairs_Bfp", Bfp, None),
    ("mol_topological_torsion_Bfp", Bfp, None),
    ("mol_pattern_Bfp", Bfp, None),
    ("mol_morgan_Bfp", Bfp, None),
    ("mol_feat_morgan_Bfp", Bfp, None),
    ("Bfp_tanimoto", Float, None),
    ("Bfp_dice", Float, None),
    ("Bfp_length", Integer, None),
    ("Bfp_weight", Integer, None),
    ("chemicalite_version", TEXT, None),
    ("rdkit_version", TEXT, None),
    ("rdkit_build", TEXT, None),
    ("boost_version", TEXT, None),
    ("rdtree_subset", BLOB, None),
    ("rdtree_tanimoto", BLOB, None),
]


# Iterate through _FUNCTIONS and create GenericFunction classes dynamically
for name, type_, doc in _FUNCTIONS:
    attributes = {"name": name, "inherit_cache": True}
    docs = []

    if doc is not None:
        docs.append(doc)

    if type_ is not None:
        attributes["type"] = type_

        type_str = f"{type_.__module__}.{type_.__name__}"
        docs.append(f"Return type: :class:`{type_str}`.")

    if len(docs) != 0:
        attributes["__doc__"] = "\n\n".join(docs)

    globals()[name] = type(name, (functions.GenericFunction,), attributes)

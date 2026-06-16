from src.drug_discovery_tools import DrugDiscoveryPlatform


platform = DrugDiscoveryPlatform(
    "data/ligand_library.csv",
    "data/target_protein.csv"
)

print("Computational Drug Discovery Platform")
print("=" * 45)

print("\nAvailable Target Proteins")
print("-" * 45)

for target in platform.get_targets():
    print("Target:", target["target_name"])
    print("Biological Role:", target["biological_role"])
    print("Disease Relevance:", target["disease"])
    print()

print("\nCandidate Ligand Ranking")
print("-" * 45)

ranked_ligands = platform.rank_ligands_by_docking_score()

for rank, ligand in enumerate(ranked_ligands, start=1):
    print(
        f"{rank}. {ligand['compound_name']} "
        f"| Docking Score: {ligand['docking_score']} kcal/mol"
    )

print("\nLipinski Rule-Based Candidates")
print("-" * 45)

lipinski_candidates = platform.filter_lipinski_candidates()

for ligand in lipinski_candidates:
    print("-", ligand["compound_name"])

print("\nSummary")
print("-" * 45)
print(platform.generate_candidate_summary())

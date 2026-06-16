import csv


class DrugDiscoveryPlatform:
    def __init__(self, ligand_file, target_file):
        self.ligand_file = ligand_file
        self.target_file = target_file
        self.ligands = self.load_csv(ligand_file)
        self.targets = self.load_csv(target_file)

    def load_csv(self, file_path):
        records = []

        with open(file_path, "r", encoding="utf-8") as file:
            reader = csv.DictReader(file)
            for row in reader:
                records.append(row)

        return records

    def get_targets(self):
        return self.targets

    def rank_ligands_by_docking_score(self):
        ranked_ligands = sorted(
            self.ligands,
            key=lambda ligand: float(ligand["docking_score"])
        )

        return ranked_ligands

    def filter_lipinski_candidates(self):
        candidates = []

        for ligand in self.ligands:
            molecular_weight = float(ligand["molecular_weight"])
            logp = float(ligand["logP"])
            donors = int(ligand["h_bond_donors"])
            acceptors = int(ligand["h_bond_acceptors"])

            if (
                molecular_weight <= 500
                and logp <= 5
                and donors <= 5
                and acceptors <= 10
            ):
                candidates.append(ligand)

        return candidates

    def generate_candidate_summary(self):
        ranked = self.rank_ligands_by_docking_score()
        lipinski_candidates = self.filter_lipinski_candidates()

        summary = "Computational Drug Discovery Summary\n"
        summary += "=" * 45 + "\n\n"

        summary += "Top Ranked Ligands by Docking Score:\n"
        for ligand in ranked:
            summary += (
                f"{ligand['compound_name']} | "
                f"Docking Score: {ligand['docking_score']} kcal/mol\n"
            )

        summary += "\nLigands Passing Basic Lipinski Criteria:\n"
        for ligand in lipinski_candidates:
            summary += f"{ligand['compound_name']}\n"

        return summary.strip()

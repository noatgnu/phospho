import re

import pandas as pd

from sequal.sequence import Sequence
from uniprot.parser import UniprotSequence

if __name__ == "__main__":
    a = pd.read_csv(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\PPM1H-PROTAC_dDIA\Results.pr_matrix.tsv", sep="\t")
    description = a[["Protein.Group", "Genes", "First.Protein.Description"]]
    d = {}
    with open(r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\fastaDatabase\Human_UniprotSP_Cano+Iso_052021.fasta", "rt") as fasta:
        current = ""
        for i in fasta:
            i = i.strip()
            if i.startswith(">"):
                u = UniprotSequence(i, True)
                if u.isotype:
                    current = f"{u.accession.replace('-', '')}{u.isotype}"
                else:
                    current = f"{u.accession.replace('-', '')}"
                d[current] = ""
            else:
                d[current] += i

    phospho_map = {}

    for i, r in a.iterrows():
        positions = []
        if r["Protein.Group"] not in phospho_map:
            phospho_map[r["Protein.Group"]] = set()
        protein = r["Protein.Group"].split(";")[0]
        seq = Sequence(r["Modified.Sequence"])

        stripped = seq.to_stripped_string()
        replaced = re.sub("^A*", "", stripped)
        ind = d[protein].index(replaced)
        a.at[i, "Peptide_start"] = ind
        a.at[i, "Peptide_end"] = ind + len(stripped)
        for n, b in enumerate(seq):
            if b.mods:
                if b.mods[0].value == "UniMod:21":

                    positions.append(str(ind+n))
                    phospho_map[r["Protein.Group"]].add(ind+n)
                    #window = d[protein][ind+n-5:ind+n+1]
        a.at[i, "Pos"] = ";".join(positions)

    result = []
    peptides = {}
    for i, r in a.iterrows():
        if r["Pos"] != "":
            positions = [int(p) for p in r["Pos"].split(";")]
        else:
            positions = []
        if phospho_map[r["Protein.Group"]]:
            if r["Protein.Group"] not in peptides:
                peptides[r["Protein.Group"]] = {}
            for s in phospho_map[r["Protein.Group"]]:
                if s not in peptides[r["Protein.Group"]]:
                    peptides[r["Protein.Group"]][s] = set()
                if s >= r["Peptide_start"] and s < r["Peptide_end"]:
                    r["Position"] = s
                    if s in positions:
                        r["Phosphorylated"] = True
                    else:
                        r["Phosphorylated"] = False
                    result.append(r.copy(deep=True))
                    peptides[r["Protein.Group"]][s].add(r["Modified.Sequence"])

    result = pd.DataFrame(result)
    columns = [i for i in result.columns if i.endswith(".raw")]
    result = result[["Protein.Group", "Position", "Phosphorylated"] + columns]
    fin = []
    for i, g in result.groupby(["Protein.Group", "Position"]):
        motif = d[i[0].split(";")[0]][i[1] - 10: i[1] + 11]
        not_phosphorylated = g[g["Phosphorylated"] == False]
        if not not_phosphorylated.empty:
            np = not_phosphorylated[columns].sum(axis=0, skipna=True)
            n = not_phosphorylated.iloc[0]
            n.loc[columns] = np
            n["Windows"] = motif
            fin.append(n)
        phosphorylated = g[g["Phosphorylated"] == True]
        if not phosphorylated.empty:
            ph = phosphorylated[columns].sum(axis=0, skipna=True)
            p = phosphorylated.iloc[0]
            p.loc[columns] = ph
            p["Windows"] = motif
            fin.append(p)

    fin = pd.DataFrame(fin)
    proportion = []
    for i, g in fin.groupby(["Protein.Group", "Position"], sort=False):
        if len(g.index) > 1:
            total = g[columns].sum(axis=0, skipna=True)
            for i2, r in g.iterrows():
                for c in columns:
                    if pd.notnull(r[c]):
                        r[c] = r[c]/total[c]
                    else:
                        r[c] = 0
                proportion.append(r)
        else:
            r = g.iloc[0]
            r[columns] = 1
            proportion.append(r)
    proportion = pd.DataFrame(proportion)
    description.drop_duplicates(inplace=True, ignore_index=True)
    for i, r in fin.iterrows():
        fin.at[i, "Modified.Sequence"] = ";".join(peptides[r["Protein.Group"]][r["Position"]])
    fin["Position"] = fin["Position"] + 1
    fin["Residue"] = fin["Windows"].str[10]
    fin = fin.merge(description, left_on="Protein.Group", right_on="Protein.Group")
    fin.to_csv("test.txt", sep="\t", index=False)
    for i, r in proportion.iterrows():
        proportion.at[i, "Modified.Sequence"] = ";".join(peptides[r["Protein.Group"]][r["Position"]])
    proportion["Position"] = proportion["Position"] + 1
    proportion["Residue"] = proportion["Windows"].str[10]
    proportion = proportion.merge(description, left_on="Protein.Group", right_on="Protein.Group")
    proportion.to_csv("test_proportion.txt", sep="\t", index=False)




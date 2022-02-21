import pandas as pd

from sequal.sequence import Sequence
from uniprot.parser import UniprotSequence, UniprotParser

if __name__ == "__main__":
    df = pd.read_csv("For_pSTY-Motif.txt", sep="\t")
    position = "N: Position"
    proteins = "T: Proteins"
    position_within = "T: Positions within proteins"
    df["protein_list"] = df[proteins].str.split(";")
    df["position_list"] = df[position_within].str.split(";")
    uni_acc = set()
    for i, r in df.iterrows():
        for p in r["protein_list"]:
            seq = UniprotSequence(p, True)
            if seq.accession:
                uni_acc.add(seq.accession.replace("-", ""))


    uniparser = UniprotParser(uni_acc, unique=True)

    with open("fasta.txt", "wt") as fasta:
        for i in uniparser.parse("fasta", "post"):
            fasta.write(i)

    d = {}

    with open("fasta.txt", "rt") as fasta:
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

    result = []
    for i, r in df.iterrows():
        for pr, po in zip(r["protein_list"], r["position_list"]):
            print(pr, po)
            if pr in d:
                a = int(po) - 1
                l = len(d[pr.strip()])

                position_in_window = a
                if a - 7 >= 0:
                    if a + 8 <= l:
                        window = d[pr][a-7:a+8]
                    else:
                        window = d[pr][a-7:l]
                    position_in_window = a - a + 7
                else:
                    if a + 8 <= l:
                        window = d[pr][0:a + 8]
                    else:
                        window = d[pr][0:l]
                r["Isoform Unique ID"] = f"{r['T: Unique identifier']};{pr};{r['T: Gene names']}"
                r["Isoform"] = pr
                r["Isoform Position"] = po
                r["Isoform Position in Window"] = position_in_window + 1
                r["Isoform Window"] = window
                result.append(r.copy())
    result = pd.DataFrame(result)

    for i, r in result.iterrows():
        p = r["Isoform Position in Window"]-1
        seq = Sequence(r["Isoform Window"])
        if p-3 >= 0:
            if seq[p-3].value in ["R", "K"]:
                if p-5 >= 0 and p+1 < len(seq):
                    if seq[p-5].value in ["V", "I", "L", "F", "W", "Y", "M"]:
                        if seq[p+1].value != "P":
                            result.at[i, "Group"] = "A"
                        else:
                            result.at[i, "Group"] = "B"
                        result.at[i, "Hydrophobic"] = seq[p - 5].value
                    else:
                        if seq[p+1].value != "P":
                            result.at[i, "Group"] = "C"
                        else:
                            result.at[i, "Group"] = "D"
    result.to_csv("test.csv", index=False)
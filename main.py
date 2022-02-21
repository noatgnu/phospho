import pandas as pd

from sequal.sequence import Sequence

if __name__ == "__main__":
    probability = "T: Phospho (STY) Probabilities"
    position = "N: Position in peptide"
    df = pd.read_csv("For_pSTY-Motif.txt", sep="\t")
    for i, r in df.iterrows():
        p = r[position]-1
        seq = Sequence(r[probability])
        if p-3 >= 0:
            if seq[p-3].value in ["R", "K"]:
                if p-5 >= 0 and p+1 < len(seq):
                    if seq[p-5].value in ["V", "I", "L", "F", "W", "Y", "M"]:
                        if seq[p+1].value != "P":
                            df.at[i, "Group"] = "A"
                        else:
                            df.at[i, "Group"] = "B"
                        df.at[i, "Hydrophobic"] = seq[p - 5].value
                    else:
                        if seq[p+1].value != "P":
                            df.at[i, "Group"] = "C"
                        else:
                            df.at[i, "Group"] = "D"

    df.to_csv("test.csv")
import re

import pandas as pd

from sequal.sequence import Sequence
from uniprot.parser import UniprotSequence

if __name__ == "__main__":
    a = pd.read_excel(r"C:\Users\toanp\Downloads\NUAK1_WT_KO_Phosphosites.xlsx")
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
    pattern = re.compile(r"[ILAVFG]\w[RK]\w\w[ST]")
    result = []
    columns = []
    for i, r in a.iterrows():
        seq = Sequence(r["Phospho (STY) Probabilities"])
        for n, a in enumerate(seq):
            if a.mods:
                probabilities = float(a.mods[0].value)
                #if probabilities >= 0.9:
                ind = d[r["Protein"]].index(seq.to_stripped_string())
                window = d[r["Protein"]][ind+n-5:ind+n+1]
                s = re.search(pattern, window)
                m = r.copy(deep=True)
                m["Probability"] = probabilities
                m["Phosphorylation Position"] = ind + n + 1
                if s:
                    m["Motif"] = window
                else:
                    m["Motif"] = ""
                result.append(m)

    result = pd.DataFrame(result)
    result.to_csv("result.txt", index=False, sep="\t")

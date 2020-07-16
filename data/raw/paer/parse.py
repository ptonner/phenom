import pandas as pd
import numpy as np
import re
import os


def parseFile(f, strain=None, header=0):

    pa = pd.read_excel(f, sheet_name=None, header=header)
    useSheets = [
        s
        for s in pa.keys()
        if re.match("(pH ?[0-9.]+,? [0-7.]+ ?(mM)?)|(All pH 0 mM)", s)
    ]

    data = meta = None

    for s in useSheets:

        m = re.match("pH ?([0-9.]+),? ([0-7.]+) ?(mM)?", s)

        if m:
            ph = m.group(1)
            la = m.group(2)
        else:
            ph = None
            la = 0

        temp = pa[s]

        if ph is None:
            temp = temp.iloc[:, :8]

            newmeta = pd.DataFrame([[ph, la]] * 7, columns=["pH", "mM-acid"])
            newmeta["pH"] = [7, 6.5, 6, 5.5, 5, 4.5, 4]
            newmeta["batch"] = 1
        else:
            temp = temp.iloc[:, :7]

            newmeta = pd.DataFrame([[ph, la]] * 6, columns=["pH", "mM-acid"])
            newmeta["batch"] = [0] * 3 + [1] * 3

        if strain is not None:
            newmeta["strain"] = strain

        if meta is None:
            meta = newmeta
        else:
            meta = pd.concat((meta, newmeta), 0)

        if data is None:
            data = temp
        else:
            if header is None:
                data = pd.merge(data, temp, on=0)
            else:
                data = pd.merge(data, temp, on="Time (hours)")

    data.columns = ["time"] + np.arange(data.shape[1] - 1).tolist()
    meta.index = range(meta.shape[0])
    data = data.set_index("time")

    return data, meta


if __name__ == "__main__":

    from collections import defaultdict

    datasets = defaultdict(list)
    for target, strain, acid, header in [
        ("PA1054 Sodium Benzoate", "PA1054", "sodium-benzoate", 0),
        ("PA1054 Citric Acid", "PA1054", "citric", 0),
        ("PA1054 Potassium Sorbate", "PA1054", "potassium-sorbate", 0),
        ("PA1054 Butyric Acid", "PA1054", "butyric", 0),
        ("PA01 Malic 09.03.17", "PA01", "malic", 0),
        ("PA01 Lactic 06.03.17", "PA01", "lactic", 0),
        ("PA01 citric 15 min time points 06.03.17", "PA01", "citric", 0),
        ("PA01 Benzoate 15 min time points", "PA01", "benzoate", 0),
        ("PA1054 Propionic Acid", "PA1054", "propionic", None),
        ("PA1054 Acetic Acid", "PA1054", "acetic", None),
        ("PA01 Potassium Sorbate 01.12.16", "PA01", "potassium-sorbate", None),
        ("PA01 Butyric Acid 15 min time points 10.11.16", "PA01", "butyric", None),
        ("Propionic acid 15 min time points PA01 02.11.16", "PA01", "propionic", 0),
        ("PA01 Acetic 15 min time points 14.10.16", "PA01", "acetic", 0),
        ("PAB Lactic Acid (1)", "PAB", "lactic", 0),
        ("PA01 Lactic Acid (1)", "PA01", "lactic", 0),
        ("PA1054 Lactic Acid", "PA1054", "lactic", 0),
        ("PA1054 Malic Acid", "PA1054", "malic", 0),
        ("PA01 Benzoate repeat 19.07.17", "PA01", "benzoate", None),
        ("PA01 Citric rerun 11.07.17", "PA01", "citric", None),
        ("PA01 Lactic repeat 13.07.17", "PA01", "lactic", None),
        ("PA01 Malic repeat 27.07.17", "PA01", "malic", None),
    ]:
        data, meta = parseFile("{}.xlsx".format(target), strain, header)

        meta["acid"] = acid
        meta["genus"] = "pseudomonas"
        meta["strain"] = strain
        meta["plate"] = target
        meta.index = range(meta.shape[0])

        d = "%s-%s" % (strain, acid)

        datasets[d].append((data, meta))

    for d, sets in datasets.items():
        data = meta = None
        for dd, m in sets:
            if data is None:
                data = dd
                meta = m
            else:
                data = pd.concat((data, dd), axis=1)
                meta = pd.concat((meta, m))

        data = data.set_axis(np.arange(data.shape[1]), axis=1)

        # remove missing data
        sel = ~data.isnull().all(0)
        data = data.loc[:, sel]
        meta = meta.loc[sel.values, :]

        assert meta.shape[0] == data.shape[1]

        print(d, data.shape, meta.shape, data.isnull().all(0).sum())

        # remove spaces from filename
        d = d.replace(" ", "_")

        path = os.path.join("../../pseudomonas", d)
        os.makedirs(path, exist_ok=True)

        data.to_csv(os.path.join(path, "data.csv"),)
        meta.to_csv(os.path.join(path, "meta.csv"), index=False)

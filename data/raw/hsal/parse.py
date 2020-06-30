import time
import datetime
import pandas as pd
from glob import glob
import os


def parse(data):

    if "Blank" in data.columns:
        del data["Blank"]

    assert "Time" in data.columns
    data["Time"] = convert_time(data["Time"])

    return data.set_index("Time")


def convert_time(time, r=2):
    time = time.apply(_parse_time)
    time = time.apply(_convert_delta_to_hours, args=(time[0],))
    time = time.round(r)
    return time


def _convert_delta_to_hours(time, t0):
    delta = time - t0
    # return delta
    return 24 * delta.days + float(delta.seconds) / 3600


def _parse_time(t):
    try:
        return datetime.datetime(*time.struct_time(time.strptime(t, "%H:%M:%S"))[:-2])
    except ValueError:
        try:
            t = time.strptime(t, "%d %H:%M:%S")
            t = list(t)
            t[2] += 1

            return datetime.datetime(*time.struct_time(t)[:-2])
        except ValueError:
            raise Exception("Time format unknown")


def save(data, key, sel, path):
    assert data.shape[1] == key.shape[0]
    data = data.loc[:, sel.values]
    key = key.loc[sel.values, :]

    data.to_csv(os.path.join(path, "data.csv"))
    key.to_csv(os.path.join(path, "meta.csv"))


if __name__ == "__main__":
    src = glob("*/")

    for d in src:
        data = glob(os.path.join(d, "*.csv"))[0]
        key = glob(os.path.join(d, "*key.xlsx"))[0]

        print("Reading in {} and {}.".format(data, key))

        data = parse(pd.read_csv(data, encoding="utf-16"))

        key = pd.read_excel(key)
        key["plate"] = d

        # print(data.head())
        # print(key.head())
        # print(data.columns)
        # print()

        # standard
        path = os.path.join("../../standard/", d)
        os.makedirs(path, exist_ok=True)
        sel = (key["mM PQ"] == 0.0) & (key.Strain == "ura3")
        if "M NaCl" in key:
            sel = sel & (key["M NaCl"] == 4.2)
        save(data, key, sel, path)

        # low
        path = os.path.join("../../low-oxidative/", d)
        os.makedirs(path, exist_ok=True)
        sel = ((key["mM PQ"] == 0.0) | (key["mM PQ"] == 0.083)) & (key.Strain == "ura3")
        if "M NaCl" in key:
            sel = sel & (key["M NaCl"] == 4.2)
        save(data, key, sel, path)

        # hi
        path = os.path.join("../../hi-oxidative/", d)
        os.makedirs(path, exist_ok=True)
        sel = ((key["mM PQ"] == 0.0) | (key["mM PQ"] == 0.333)) & (key.Strain == "ura3")
        if "M NaCl" in key:
            sel = sel & (key["M NaCl"] == 4.2)
        save(data, key, sel, path)

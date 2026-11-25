"""Normalize the NIST StRD linear-regression files (tassel-6#183, ADR-0029).

Usage: nist.py IN_DIR OUT_DIR

For each dataset writes OUT_DIR/<Name>.data.tsv (header row, then the data rows with the
numbers exactly as NIST prints them) and OUT_DIR/<Name>.certified.txt (`key = value` lines,
every certified number copied as NIST prints it, never rounded or reformatted). Runs in the
pinned Python container; standard library only.
"""
import re
import sys
from pathlib import Path

DATASETS = [
    "Norris", "Pontius", "NoInt1", "NoInt2", "Filip", "Longley",
    "Wampler1", "Wampler2", "Wampler3", "Wampler4", "Wampler5",
]
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?"


def parse(name, text):
    lines = text.splitlines()
    data_at = next(i for i, l in enumerate(lines) if re.match(r"Data:\s+y(\s|$)", l))
    head, body = lines[:data_at], lines[data_at:]
    certified = {"dataset": name}

    observations = re.search(r"(\d+) Observations", "\n".join(head))
    certified["observations"] = observations.group(1)
    model_at = next(i for i, l in enumerate(head) if l.startswith("Model:"))
    model_lines = []
    for l in head[model_at + 1:]:
        if "y =" in l or model_lines:
            if not l.strip():
                break
            model_lines.append(l.strip())
    model = " ".join(model_lines)
    certified["model"] = re.sub(r"^y\s*=\s*", "", re.sub(r"\s*\+\s*e$", "", model))
    certified["intercept"] = "false" if name.startswith("NoInt") else "true"

    parameters = []
    for l in head:
        m = re.match(rf"\s*(B\d+)\s+({NUM})\s+({NUM})\s*$", l)
        if m:
            parameters.append(m.groups())
    certified["parameters"] = str(len(parameters))
    for key, estimate, sd in parameters:
        certified[key] = f"{estimate} {sd}"

    for l in head:
        m = re.match(rf"\s*Standard Deviation\s+({NUM})\s*$", l)
        if m:
            certified["residual_sd"] = m.group(1)
        m = re.match(rf"\s*R-Squared\s+({NUM})\s*$", l)
        if m:
            certified["r_squared"] = m.group(1)
        m = re.match(rf"Regression\s+(\d+)\s+({NUM})\s+({NUM})\s+({NUM}|Infinity)\s*$", l)
        if m:
            (certified["df_regression"], certified["ss_regression"],
             certified["ms_regression"], certified["f_statistic"]) = m.groups()
        m = re.match(rf"Residual\s+(\d+)\s+({NUM})\s+({NUM})\s*$", l)
        if m:
            certified["df_residual"], certified["ss_residual"], certified["ms_residual"] = m.groups()

    names = body[0].split()[1:]
    rows = [l.split() for l in body[1:] if l.strip()]
    assert all(len(r) == len(names) for r in rows), name
    assert len(rows) == int(certified["observations"]), name
    for expected in ("residual_sd", "r_squared", "df_regression", "df_residual"):
        assert expected in certified, (name, expected)
    return certified, names, rows


def main(in_dir, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name in DATASETS:
        text = (Path(in_dir) / f"{name}.dat").read_text(encoding="ascii")
        certified, names, rows = parse(name, text)
        (out / f"{name}.certified.txt").write_text(
            "".join(f"{k} = {v}\n" for k, v in certified.items()), encoding="ascii")
        (out / f"{name}.data.tsv").write_text(
            "\t".join(names) + "\n" + "".join("\t".join(r) + "\n" for r in rows), encoding="ascii")


if __name__ == "__main__":
    main(*sys.argv[1:3])

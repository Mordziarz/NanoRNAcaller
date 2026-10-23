#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_samples_table(filename):
    with open(filename, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        required = {"sample", "group"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "samples table must contain sample and group columns"
            )

        rows = list(reader)

    groups = {}

    for row in rows:
        sample = (row.get("sample") or "").strip()
        group = (row.get("group") or "").strip()

        if not sample or not group:
            raise ValueError(
                "Each sample must have a non-empty sample and group"
            )

        groups.setdefault(group, []).append(sample)

    if len(groups) != 2:
        raise ValueError(
            "SUPPA diffSplice requires exactly two groups here. "
            "Found groups: " + ", ".join(groups.keys())
        )

    return groups


def split_matrix(filename, groups, output_dir, suffix):
    with open(filename, encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)

        if not header or len(header) < 1:
            raise ValueError(
                "File {} does not contain a valid header".format(filename)
            )

        header = [x.strip() for x in header]

        all_known_samples = [s for samples in groups.values() for s in samples]
        if header[0] in all_known_samples:
            header = [""] + header

        sample_columns = {}

        for group_name, samples in groups.items():
            missing = [s for s in samples if s not in header[1:]]

            if missing:
                raise ValueError(
                    "File {} is missing samples {}. "
                    "Available columns: {}".format(
                        filename, missing, header[1:]
                    )
                )

            sample_columns[group_name] = [
                header.index(sample) for sample in samples
            ]

        output_handles = {}
        output_writers = {}

        try:
            for index, group_name in enumerate(groups.keys(), start=1):
                output_path = output_dir / (
                    "group{}.{}".format(index, suffix)
                )

                output_handles[group_name] = open(
                    output_path,
                    "w",
                    encoding="utf-8",
                    newline=""
                )

                writer = csv.writer(
                    output_handles[group_name],
                    delimiter="\t",
                    lineterminator="\n"
                )

                columns = sample_columns[group_name]

                writer.writerow(
                    [header[0]] + [header[i] for i in columns]
                )

                output_writers[group_name] = writer

            for row in reader:
                if not row:
                    continue

                for group_name in groups:
                    columns = sample_columns[group_name]

                    values = [
                        row[i] if i < len(row) else "" 
                        for i in columns
                    ]

                    output_writers[group_name].writerow(
                        [row[0]] + values
                    )

        finally:
            for handle in output_handles.values():
                handle.close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--psi", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--outdir", required=True)

    args = parser.parse_args()

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    groups = read_samples_table(args.samples)

    group_names = list(groups.keys())

    print("Compared groups:")
    print("  group1 =", group_names[0], groups[group_names[0]])
    print("  group2 =", group_names[1], groups[group_names[1]])

    split_matrix(
        args.psi,
        groups,
        output_dir,
        "psi"
    )

    split_matrix(
        args.expression,
        groups,
        output_dir,
        "tpm"
    )

    print("Group files saved in:", output_dir)


if __name__ == "__main__":
    main()

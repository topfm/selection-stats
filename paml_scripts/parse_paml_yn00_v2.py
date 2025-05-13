import csv

def write_results(results, out_file):
    rows = []
    for iso in results:
        for comp_iso in results[iso]:
            row = {
                "isolate": iso,
                "comparedTo_iso": comp_iso
            }
            for method in results[iso][comp_iso]:
                for stat_key, stat_val in results[iso][comp_iso][method].items():
                    col_name = f"{method}_{stat_key}"
                    row[col_name] = stat_val
            rows.append(row)

    # Write to file with consistent headers
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=sorted(set().union(*[r.keys() for r in rows])), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

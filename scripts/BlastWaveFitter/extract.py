import ROOT
import csv
import os

# Load mapping from CSV
csv_path = "graph_mapping.csv"  # Adjust path as needed
output = ROOT.TFile("spec_lambda_infered.root", "RECREATE")

with open(csv_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    current_file = None
    current_root = None
    for row in reader:
        file = row["file"]
        table = row["table"]
        graph_name = row["graph"]
        cent = row["cent"]
        category = row["category"]

        if current_file != file:
            if current_root:
                current_root.Close()
            print(f"Opening {file}...")
            current_root = ROOT.TFile.Open(file)
            if not current_root or current_root.IsZombie():
                print(f"Failed to open {file}")
                continue
            current_file = file

        dir = current_root.GetDirectory(table)
        if not dir:
            print(f"Warning: {table} not found in {file}")
            continue
        graph = dir.Get(graph_name)
        if not graph or not isinstance(graph, (ROOT.TGraphErrors, ROOT.TGraphAsymmErrors)):
            print(f"Warning: {graph_name} not found or not TGraphErrors or TGraphAsymmErrors in {file}/{table}")
            continue

        graph.SetName(f"{category}_{cent}")
        output.cd()
        graph.Write()
        print(f"Wrote {category}_{cent}")

    if current_root:
        current_root.Close()

output.Close()
print("Done. Extracted graphs are in spec_lambda_infered.root")

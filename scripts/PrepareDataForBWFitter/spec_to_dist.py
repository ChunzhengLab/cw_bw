


#!/usr/bin/env python
"""
Script: spec_to_dist.py

Purpose:
    Convert TGraphAsymmErrors objects named "spec_*" in an input ROOT file
    into TH1D histograms representing pT distributions (value = x * y) with
    variable bin edges aligned to graph points.

Usage:
    python spec_to_dist.py input_file.root

Input:
    input_file.root - ROOT file containing TGraphAsymmErrors named "spec_*".

Output:
    <basename>_spec_to_dist.root - ROOT file with TH1D histograms named
    "pT_dist_<suffix>".
"""

import sys
import os
from array import array
import ROOT

def main():
    if len(sys.argv) != 2:
        print("Usage: {} input_file.root".format(sys.argv[0]))
        sys.exit(1)
    input_file = sys.argv[1]
    basename = os.path.splitext(os.path.basename(input_file))[0]
    outfile_name = "{}_spec_to_dist.root".format(basename)

    # Open input and output ROOT files
    fin = ROOT.TFile.Open(input_file, "READ")
    if not fin or fin.IsZombie():
        print("Error: cannot open input file", input_file)
        sys.exit(1)
    fout = ROOT.TFile.Open(outfile_name, "RECREATE")
    
    # Loop over all objects in the file
    for key in fin.GetListOfKeys():
        obj = key.ReadObj()
        # Select TGraphAsymmErrors objects whose name starts with "spec_"
        if obj.InheritsFrom("TGraphAsymmErrors") and obj.GetName().startswith("spec_"):
            name = obj.GetName()             # e.g. "spec_xxx"
            suffix = name[len("spec_"):]     # "xxx"
            n = obj.GetN()

            # Retrieve x and y arrays
            x_arr = obj.GetX()
            y_arr = obj.GetY()

            # Build variable bin edges so histogram bins align with x points
            if n > 1:
                mid_edges = [(x_arr[i] + x_arr[i+1]) / 2.0 for i in range(n-1)]
                low_edge  = x_arr[0] - (mid_edges[0] - x_arr[0])
                high_edge = x_arr[-1] + (x_arr[-1] - mid_edges[-1])
                bin_edges = array('d', [low_edge] + mid_edges + [high_edge])
            else:
                # single-point graph: give a small default width
                width = x_arr[0] * 0.1 if x_arr[0] != 0 else 1.0
                bin_edges = array('d', [x_arr[0] - width/2, x_arr[0] + width/2])

            # Create TH1D and fill with x*y values
            hist_name = "pT_dist_" + suffix
            hist = ROOT.TH1D(hist_name, "pT distribution from {}".format(name), n, bin_edges)
            for i in range(n):
                value = x_arr[i] * y_arr[i]
                hist.SetBinContent(i+1, value)

            # Write histogram to output file
            fout.WriteTObject(hist, hist.GetName())

    fout.Close()
    print("Wrote output file:", outfile_name)

if __name__ == "__main__":
    main()
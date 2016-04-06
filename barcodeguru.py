#!/usr/bin/env python

# barcodeguru.py
#
# Copyright 2013-2015 James S Blachly, MD and The Ohio State University
#
# License: licensed under GNU GPL 3 ; see LICENSE

import sys
import re
import numpy as np

import rules

def main():
    print("Enter barcode sequences (e.g. AACCTG; do not use quotes) in one of the following ways:")
    print("\t1) one per line")
    print("\t2) as a single line comma-separated list (NOT IMPLEMENTED)")
    print("Press CTRL-D twice when finished.\n(Windows users try CTRL-Z)\n")

    barcodes = []
    for line in sys.stdin:
        line = line.upper().strip('\r\n')
        barcodes.append(line)
    
    print("")
    if not rules.checkUniformLength(barcodes):
        print("Index lengths:\t FAILED. Index length mismatch(es) detected.")
    print("Duplication:\t" + rules.Duplication(barcodes))
    print("Pool size:\t" + rules.PoolSize(barcodes))
    print("Dual Lasers:\t" + rules.Lasers(barcodes))
    print("Base Matches:\t" + rules.BaseMatch(barcodes))
    print("Nucleotide frequencies:\t" + rules.pwm(barcodes))
    print("")
    return True


if __name__ == "__main__":
    main()


#!/usr/bin/python

# rules.py
#
# Copyright 2013-2015 James S Blachly, MD and The Ohio State University
#
# License: licensed under GNU GPL 3 ; see LICENSE

import sys
import re
import numpy as np

def checkUniformLength(barcodes):
    """Check that all indexes have same length.
    Used by rule functions.
    Returns Bool."""
    
    barcodeLenPrior = -1
    
    # ensure all indexes have same length
    for barcode in barcodes:
        if barcodeLenPrior == -1:
            barcodeLenPrior = len(barcode)
        else:
            if len(barcode) != barcodeLenPrior:
                return False
    return True

def Duplication(barcodes):
    """Ensure that each index appears only once.

    returns string."""
    
    if len(barcodes) == len(set(barcodes)):
        return "OK"
    else:
        return "FAILED: Index duplication detected."

def PoolSize(barcodes):
    """Number of indexes must be 1 (unpooled) or greater than 3 (ie 4 or more pooled samples)

    returns string."""
    
    if len(barcodes) == 1:
        return "OK"
    elif len(barcodes) > 3:
        return "OK"
    else:
        return "FAILED: no. of samples/lane cannot be 2 or 3"

def Lasers(bc):
    """Ensure that for each base position, both lasers are firing.

    returns string."""
    
    # this function modifies barcodes with a regex
    # A/C->M
    # G/T->K
    # so I must make a copy of barcodes so as not to alter the original

    barcodes = list(bc)

    if not checkUniformLength(barcodes):
        return "FAILED: Cannot test until index length mismatch(es) fixed."

    # I would have said for index in barcodes:,
    # but this gets confusing as I am also using
    # Python's list.index() function ...
    for barcode in barcodes:
        n = barcodes.index(barcode)
        barcode = re.sub('[AC]', 'M', barcode)  # IUPAC M = A or C
        barcode = re.sub('[GT]', 'K', barcode)  # IUPAC K = G or T
        barcodes[n] = barcode
    
    # the check for uniform index length necessary prior to this step
    # as it presupposes identical index lengths
    for basePos in range(len(barcodes[0])):
        bases = []
        for barcode in barcodes:
            bases.append(barcode[basePos])
        if len(set(bases)) == 1:        # if both M + K at position set cardinality is 2
            return "FAILED: Each base position must have both {A or C} and {G or T}."
    return "OK"

def BaseMatch(barcodes):
    """
    Ensure appropriate Hamming distance between any two barcodes
    Appropriate distance is > 2 to prevent barcode collisions 
    when allowing 1 mismatch during demuxing
    Example:
    Consider the barcodes AAAGGG / ATTGGG. They have a distance/similarity
    of 2. If the instrument reads AATGGG, it is equidistant from both of them

    AAAGGG <- AATGGG -> ATTGGG

    For standard 6 nt barcodes, this corresponds to requiring a similarity < 4
    (i.e. erroring if similarity > 3)
    
    Because Hamming distance = length - similarity, we will measure similarity
    and return "OK" if distance > 2. 

    returns string."""

    if not checkUniformLength(barcodes):
        return"FAILED: Cannot test until index length mismatch(es) fixed."
    
    numBC = len(barcodes)
    simMatrix = np.asmatrix(np.zeros((numBC,numBC), np.int32))
    
    # loop through the upper triangular matrix
    for i in range(numBC):
        for j in range(i,numBC):
            if j == i:
                pass
            else:
                barcode1 = barcodes[i]
                barcode2 = barcodes[j]
                simMatrix[i,j] = similarity(barcode1, barcode2)
                simMatrix[j,i] = simMatrix[i,j] # make symmetric
    
#   for barcode1 in barcodes:
#       i = barcodes.index(barcode1)
#       for barcode2 in barcodes:
#           j = barcodes.index(barcode2)
#           if j == i:
#               pass
#           else:
#               simMatrix[i,j] = similarity(barcode1, barcode2)
    
    #print(simMatrix)
    
    if (len(barcodes[0]) - simMatrix.max()) <= 2:
        failStr = "FAILED: The following indexes were too similar ( Hamming distance <= 2 ):\n"
        failStr += "\t\tBC1\tBC2\tSimilarity\n"
        for i in range(numBC):
            for j in range(i,numBC):
                if simMatrix[i,j] > 3:
                    failStr += ('\t\t' + barcodes[i] + '\t' + barcodes[j] + '\t' + str(simMatrix[i,j]) + '\n' )
        return failStr
    else:
        return "OK"

def similarity(barcode1, barcode2):
    """Compute a measure of similarity between barcode1 and barcode 2:
    
    Counts the number of times barcode1 and barcode 2 have the same base
    at position i. 
    (length - similarity = Hamming distance)

    returns an integer."""
    
    simscore = 0
    for i in range(len(barcode1)):      # NB req uniform length
        if barcode1[i] == barcode2[i]:
            simscore += 1
    
    return simscore
    
def entropy(barcodes):
    """Compute Shannon entropy (in bits) of supplied barcode

    returns an integer
    """
    if not checkUniformLength(barcodes):
        return"FAILED: Cannot test until index length mismatch(es) fixed."
    
    numBC = len(barcodes)
    simMatrix = np.asmatrix(np.zeros((numBC,numBC), np.int32))
    return True

def pwm(barcodes):
    """Print the position weight matrix
    
    returns string"""

    if not checkUniformLength(barcodes):
        return"FAILED: Cannot test until index length mismatch(es) fixed."
    
    numBCs = len(barcodes)
    lenBCs = len(barcodes[0])
    PFMatrix = np.asmatrix(np.zeros((4,lenBCs), np.float))

    for i,nt in enumerate(['A','C','G','T']):
        for j in xrange(lenBCs):
            PFMatrix[i,j] = count_nts_at_pos(barcodes, nt, j)

    print "Position weight matrix: (A/C/G/T)"
    PWMatrix = PFMatrix / np.float(numBCs)

    print PWMatrix
    if PWMatrix.min() <= 0.1:
        return "WARNING: entropy is low. At least one position has a nt freq <= 10%"
    else:
        return "OK"

def count_nts_at_pos(barcodes, nt, pos):
    sum = 0    
    for bc in barcodes:
        if bc[pos] == nt: sum = sum+1
    return sum

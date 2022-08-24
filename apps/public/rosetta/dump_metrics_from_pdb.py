#!/usr/bin/env python3

#Author Jared Adolf-Bryfogle
#Seems to be broken.

from argparse import ArgumentParser
from jade2.pymol_jade.PyMolScriptWriter import *
from jade2.rosetta_jade.ScoreFiles import ScoreFile
from collections import defaultdict
from jade2.basic.figure.creation import *
import pandas
import os, ast


skip_types=[
    "ATOM",
    "HETATOM",
    "REMARK",
    "TER",
    "HELIX",
    "SHEET",
    "EXPDTA",
    "HEADER",
    "SSBOND",
    "LINK",
    "SEQRES",
    "AUTHOR",
    "KEYWDS",
    "SOURCE",
    "CAVEAT",
    "COMPND",
    "MODEL",
    "LINK",
    "CONECT"
]

def deduce_str_type(s):
    """
    Deduce the type of a string.  Either return the string as the literal, or as the string if not possible.
    http://stackoverflow.com/questions/13582142/deduce-the-type-of-data-in-a-string

    :param s: str
    :return:
    """
    try:
        return ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return s

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Extracts metric data from a list of PDBs and writes a scorefile in json format.")

    parser.add_argument("-l", "--pdblist",
                        help = "List of PDBs to use",
                        default = "PDBLIST.txt")

    parser.add_argument("-o", "--scorefile",
                        help = "Output scorefile name",
                        default = "score.json")

    options = parser.parse_args()

    pdbfiles = []
    for line in open(options.pdblist, 'r'):
        line = line.strip()
        if not line or line.startswith("#"): continue
        pdbfiles.append(line)

    all_data = []
    datatypes = []
    for pdbfile in pdbfiles:
        d = defaultdict()
        read = True
        for line in open(pdbfile, 'r'):
            line = line.strip()
            if not line: continue
            lineSP = line.split()
            if lineSP[0] in skip_types:
                continue
            
            if lineSP[0] == "#BEGIN_POSE_ENERGIES_TABLE":
                read = False

            if lineSP[0] == "#END_POSE_ENERGIES_TABLE":
                read = True

            if not read: continue
            datatype = lineSP[0]
            if len(lineSP) == 2:
                data = lineSP[1]
            else:
                data = " ".join(lineSP[0:])

            if data == "nan":
                data = "NaN"

            d[ datatype ] = deduce_str_type(data)
            if not datatype in datatypes:
                datatypes.append(datatype)

        d["decoy"] = pdbfile
        d["total_score"] = 0
        all_data.append(d)

    df = pandas.DataFrame.from_dict(all_data)
    df = detect_numeric(df)
    #df.to_csv("debugging.csv", sep=",")

    df = df.sort_values("decoy")
    df.to_json(options.scorefile, orient='records', lines=True)

    print("Datatypes: "+repr(datatypes))
    print("\n\n")
    print("Finished Writing " + options.scorefile)

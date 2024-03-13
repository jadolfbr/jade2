#!/usr/bin/env python3


import os,json,re,glob,sys
from argparse import ArgumentParser
#from jade2.rosetta_jade.score_util import parse_decoy_scores
from collections import defaultdict
import gzip


def get_pdbs(argu):
    if os.path.isdir(argu):
        print("Gathering PDBs: " + argu)
        pdbs = glob.glob(argu+"/*.pdb*")
        return pdbs
    elif os.path.isfile(argu) and not re.search(".pdb", argu) and not re.search(".pdb.gz", argu):
        print("Parsing PDBs: " + argu)
        return [ x.strip() for x in open(argu, 'r').readlines() if not x.startswith('#') and x ]
    else:
        return [argu]


def open_file(file_path, mode='r'):
    """
    Open a file which can be gzipped.

    :param file_path:
    :param mode:
    :return:
    """
    if file_path.split(".")[-1] == "gz":
        # print "opening gzipped file"
        INFILE = gzip.open(file_path, mode + 'b')
    else:
        INFILE = open(file_path, mode)

    return INFILE

def get_decoy_name(decoy):
    """
    Get the decoy name from path or name, whether .pdb, .pdb.gz or no extension.
    :param decoy:
    :rtype:str
    """

    name = os.path.basename(decoy)

    if re.search(".pdb.gz", name) or re.search(".cif.gz", name):
        return '.'.join(name.split(".")[0:-2])
    elif re.search(".pdb", name) or re.search(".cif", name):
        return '.'.join(name.split('.')[0:-1])
    else:
        return name

def parse_decoy_scores(decoy_path):
    """
    Parse a score from a decoy and return a dictionary.
    :param decoy_path:
    :return: defaultdict
    """

    data = defaultdict()
    labels = []
    scores=[]
    INFILE = open_file(decoy_path)
    for line in INFILE:
        if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
            lineSP = line.split()
            #if len(lineSP) == 1:
            data['decoy'] = get_decoy_name( os.path.basename(decoy_path) )
            #else:
            #    data['decoy'] = get_decoy_name( lineSP[-1] )

        elif line.startswith("label"):
            labels = line.split()[1:]
            labels = ["total_score" if x=="total" else x for x in labels ]

        elif line.startswith("pose"):
            scores = [ float(x) for x in line.split()[1:] ]

        else:
            continue

    INFILE.close()

    for i in range(0, len(labels)):
        data[labels[i]] = scores[i]

    return data

def get_parser():
    parser = ArgumentParser(description="This script creates a Rosetta score file from a set of structures - by parsing the score from them. Pass a directory, a PDBLIST, and/or a list of filenames")


    parser.add_argument("--prefix",
                        help = "Any prefix to use.  ",
                        default = "")

    parser.add_argument("decoys",
                        help = "A directory, a PDBLIST, and/or a list of filenames",
                        default = [],
                        nargs="*")
    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()

    #print(options)

    if len(options.decoys) == 0:
        sys.exit("Please pass decoys to parse score.")

    decoys = []
    for argu in options.decoys:

        pdbs = get_pdbs(argu)
        #print(pdbs)
        decoys.extend(pdbs)

    #print("\n".join(decoys))

    if options.prefix:
        OUTFILE = open(options.prefix+"score.json", 'w')
    else:
        OUTFILE = open(options.prefix + "score.json", 'w')

    scores = []
    decoy_num = 1
    print("Reading",len(decoys), "decoys")
    for decoy in decoys:
        if decoy_num % 50 == 0:
            print("Decoy",decoy_num)
        score_dict = parse_decoy_scores(decoy)
        if not score_dict:
            print("decoy", decoy, "has no score")
        if score_dict:
            OUTFILE.write(json.dumps(score_dict, sort_keys=True)+"\n")

        decoy_num+=1

    #OUTFILE.write("\n".join(json.dumps(scores)))

    OUTFILE.close()
    print("Done")
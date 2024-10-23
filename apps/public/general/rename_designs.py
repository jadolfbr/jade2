#!/usr/bin/env python3


import re
import sys
import os

from argparse import ArgumentParser

from jade2.basic import path

def get_parser():
    parser = ArgumentParser(description="Renames original files to new names for design ordering.  Copy all models going to be ordered into a single directory first. Run from directory with pdb files already copied in.  First, use the -n option and pass a pdblist For example -n AM will then have new names be from AM1-X.  This will create a file that will have new names to old names.  After that, run with the -i option.  This is a two-step process to make sure the names are correct.")
    parser.add_argument("-i", "--new_names",
                        help = "File with new to old names.  Example line: new_name  *  filename.  Can have lines that don't have all three.  Will only rename if it has a star in the second column.",
                        required = False)

    parser.add_argument("-p", "--prefix",
                        default="PDBLIST.txt",
                        help = "Add a prefix to input files. Must be run from same directory.")

    parser.add_argument("-l", "--pdblist",
                        help = "Files to add a prefix to. Must be run from same directory.")

    parser.add_argument('--output_numbered_name','-n',
                        help = "Write a new names txt file that is ordered from 1-N with the given new name.")

    return parser



if __name__ == "__main__":

    parser = get_parser()

    options = parser.parse_args()

    if not options.new_names and not options.output_numbered_name:
        if options.prefix and options.pdblist:
            pass
        else:
            sys.exit("If no new name file, prefix and pdblist must be passed")

        for line in open(options.pdblist, 'r').readlines():
            if not line: continue
            if line.startswith("#"): continue
            line = line.strip()
            new_name = options.prefix+line
            os.system('cp '+line+' '+new_name)
    elif options.output_numbered_name:
        if not options.pdblist: sys.exit("Must give PDB list for output as numbered name!")
        OUT = open('new_names.txt', 'w')
        old_names = open(options.pdblist, 'r').readlines()
        i = 0
        for old_name in old_names:
            old_name = old_name.strip()
            if not old_name: continue
            if old_name.startswith('#'): continue

            i+=1
            new_name = options.output_numbered_name+str(i).zfill(3)
            #print(new_name)
            OUT.write(new_name+' * '+old_name+"\n")
        OUT.close()
    else:
        INFILE = open(options.new_names, 'r')
        for line in INFILE:
            line = line.strip()
            lineSP = line.split()
            if not line or len(lineSP) != 3:
                continue

            new_name = lineSP[0].strip()
            create = lineSP[1].strip()
            old_name = lineSP[2].strip()

            if create != '*': continue

            print("Renaming "+old_name+" to "+new_name+path.get_decoy_extension(old_name))

            os.system("cp "+old_name+" "+new_name+path.get_decoy_extension(old_name))

        INFILE.close()




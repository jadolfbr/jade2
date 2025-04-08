#!/usr/bin/env python3
# coding: utf-8

#from sequence import fasta
import argparse
import re,sys
from collections import defaultdict

from jade2.basic import path
from jade2.basic.structure import util
from jade2.basic import general
from jade2.basic.structure.Structure import AntibodyStructure
from jade2.basic.structure.BioPose import BioPose
import pandas

gfp = """
SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL
VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV
NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLAD
HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
"""


igg_types = ["IgG_order", "IgG_order_lambda", "IgG_order_kappa", "IgG_order_heavy"]
format_types = ["basic", "fasta", "general_order", "csv"]
#format_types.extend(igg_types)


#Biopython removed that function. So screw biopython.
def three_to_one(three_letter_code):
    amino_acid_dict = {
        'Ala': 'A',
        'Arg': 'R',
        'Asn': 'N',
        'Asp': 'D',
        'Cys': 'C',
        'Gln': 'Q',
        'Glu': 'E',
        'Gly': 'G',
        'His': 'H',
        'Ile': 'I',
        'Leu': 'L',
        'Lys': 'K',
        'Met': 'M',
        'Phe': 'F',
        'Pro': 'P',
        'Ser': 'S',
        'Thr': 'T',
        'Trp': 'W',
        'Tyr': 'Y',
        'Val': 'V'
    }
    return amino_acid_dict.get(three_letter_code, None)

def one_to_three(one_letter_code):
    amino_acid_dict = {
        'A': 'Ala',
        'R': 'Arg',
        'N': 'Asn',
        'D': 'Asp',
        'C': 'Cys',
        'Q': 'Gln',
        'E': 'Glu',
        'G': 'Gly',
        'H': 'His',
        'I': 'Ile',
        'L': 'Leu',
        'K': 'Lys',
        'M': 'Met',
        'F': 'Phe',
        'P': 'Pro',
        'S': 'Ser',
        'T': 'Thr',
        'W': 'Trp',
        'Y': 'Tyr',
        'V': 'Val'
    }
    return amino_acid_dict.get(one_letter_code, None)


def get_parser():

    parser = argparse.ArgumentParser(description="Uses Biopython to print sequence information.  Example:\n"
                                                 "get_seq.py --pdb 2j88_A.pdb --format fasta --outpath test.txt")

    parser.add_argument("--pdb", "-s",
                        help = "Input PDB path")

    parser.add_argument("--pdblist", "-l",
                        help = "Input PDB List")

    parser.add_argument("--pdblist_input_dir", "-i",
                        help = "Input directory if needed for PDB list")

    parser.add_argument("--chain", "-c",
                        help = "A specific chain to output",
                        default = "")

    parser.add_argument("--cdr",
                        help = "Pass a specific CDR to output alignments of.",
                        default = "")

    parser.add_argument("--format",
                        help = "The output format requried.",
                        default = "fasta",
                        choices = format_types)

    parser.add_argument("--outpath", "-o",
                        default="sequences.fasta",
                        help = "Output path.  ")

    parser.add_argument("--prefix", "-t",
                        default = "",
                        help = "Tag to add before chain")

    parser.add_argument("--region",
                        help = "specify a particular region, start:end:chain")

    #parser.add_argument("--strip_n_term",
    #                    help = "Strip this sequence off the N-term of resulting sequences. (Useful for antibodies")

    parser.add_argument("--strip_c_term",
                        help = "Strip this sequence off the C-term of resulting sequences. (Useful for antibodies")

    parser.add_argument("--pad_c_term",
                        help = "Pad this sequence with some C-term (Useful for antibodies")

    parser.add_argument("--pad_to_n_residues",
                        help = "If padding c term, pad up to this number of residues.")

    parser.add_argument("--pad_seq_with_stop_codons",
                        help = "Pad N and C of the padding with stop codons",
                        default = False,
                        action="store_true")

    parser.add_argument("--output_original_seq",
                        default = False,
                        action = "store_true",
                        help = "Output the original sequence and the striped seqeunce if stripped.  Default FALSE.")

    parser.add_argument("--correct_for_tev", '-v',
                        default = False,
                        action='store_true',
                        help = "Correct for Proline TEV protease cleavage site at second position of N-terminal, by adding a G as first residue.")

    parser.add_argument("--skip_chain_output",
                        default = False,
                        action = "store_true",
                        help = "Skip outputting the chain in the fasta.  Useful for designs.")

    parser.add_argument("--three_letter",
                        default = False,
                        action = "store_true",
                        help = "Write 3 letter codes instead of one.")
    return parser

if __name__ == "__main__":



    parser = get_parser()
    options = parser.parse_args()

    if options.format in igg_types:
        igg_type_format = True
    else:
        igg_type_format = False


    options.format = options.format.lower()

    #if not options.pdb:
    #    options.pdb = sys.argv[1]

    if (options.format == "igg_order_lambda" or options.format == "igg_order_kappa") and not options.chain:
        options.chain = "L"

    if options.format == "igg_order_heavy" and not options.chain:
        options.chain = "H"


    sequences = defaultdict()
    ordered_ids = []

    pdbs = []
    if options.pdb:
        pdbs.append(options.pdb)

    if options.pdblist:
        INFILE = open(options.pdblist, 'r')
        for line in INFILE:
            line = line.strip()
            if not line or re.search("PDBLIST", line) or line[0]=='#':
                continue

            if options.pdblist_input_dir:
                pdb_path = options.pdblist_input_dir+"/"+line
            else:
                pdb_path = line
            pdbs.append(pdb_path)
        INFILE.close()

    csv_data = []
    for pdb in pdbs:
        #print "Reading "+pdb
        biopose = BioPose(pdb)
        biostructure = biopose.structure()
        ab_info = AntibodyStructure()
        if options.cdr:
            seq = ab_info.get_cdr_seq(biopose, options.cdr.upper())

            sequences[path.get_decoy_name(pdb)+"_"+options.cdr] = seq
            ordered_ids.append(path.get_decoy_name(pdb)+"_"+options.cdr)
        else:
            d = dict()
            d['path'] = pdb
            d['ID'] = path.get_decoy_name(pdb)
            for biochain in biostructure[0]:
                if options.chain and biochain.id != options.chain:
                    continue
                if util.get_chain_length(biochain) == 0:
                    continue
                seq = util.get_seq_from_biochain(biochain)
                if not seq:
                    print("Sequence not found for: ", path.get_decoy_name(pdb))
                    continue

                if  options.correct_for_tev and seq[1] == 'P':
                    seq = 'G'+seq

                #b = os.path.basename(pdb).replace('.pdb.gz', '')
                out_id = ""
                if options.skip_chain_output:
                    out_id = path.get_decoy_name(pdb)
                else:
                    out_id = path.get_decoy_name(pdb)+"_"+biochain.id
                sequences[out_id] = seq
                ordered_ids.append(out_id)
                d[biochain.id] = seq
            csv_data.append(d)




    outlines = []
    i = 1

    if options.format == "general_order" or igg_type_format:
        #outlines.append(begin_schief_order)

        print("DOUBLE CHECK FOR MISSING DENSITY IN THE STRUCTURE!  MAKE SURE SEQUENCES MATCH!")

    if options.format == "igg_order_heavy":
        outlines.append("###################################################################################################################################################")
        outlines.append("Use pFUSEss-CHIg-hG1 (human heavy) vector and clone between the two flanking regions: "
                        "GCACTTGTCACGAATTCG--[Antibody-Heavy-FV]--GCTAGCACCAAGGGCCCATC")
        outlines.append("###################################################################################################################################################\n")

        print("Sequence should flank: ANY - WGqgtlVtvss\n")
        print(">Heavy example:\nEVQLVESGGGLVKPGGSLRLSCSASGFDFDNAAMSWVRQPPGKGLEWVGTITGPGEGWSVDYAAPVEGRFTISRLNSINFLYLEMNNLRMEDSGLYFCARGEWEFRNGETYSALTYWGRGTLVTVSS")

    elif options.format == "igg_order_lambda":
        outlines.append("#################################################################################################################################################")
        outlines.append("Use pFUSEss-CLIg-hL2 (human lambda) vector and clone between the two flanking regions: "
                        "GCACTTGTCACGAATTCG [Antibody-Lambda-Fv] CTAGGTCAGCCCAAGGCT")
        outlines.append("#################################################################################################################################################\n")

        print("Sequence should flank: ANY - gsgtqvTV\n")
        print(">Lambda example:\nSYELTQETGVSVALGRTVTITCRGDSLRYHYASWYQKKPGQAPILLFYGKNNRPSGVPDRFSGSASGNRASLTISGAQAEDDAEYYCMSAAKPGSWTRTFGGGTKLTVL")

    elif options.format == "igg_order_kappa":
        outlines.append("##########################################################################################################################################################")
        outlines.append("Use pFUSEss-CLIg-hk (human kappa) vector and clone between the two flanking regions: "
                        "GCACTTGTCACGAATTCA--[Antibody-Kappa-Fv]--CGTACGGTGGCTGCACCA")
        outlines.append("##########################################################################################################################################################\n")

        print("Sequence should flank: ANY - FGGGtkveik\n")
        print(">Kappa example:\nDIQMTQSPASLSASVGETVTITCRASENIYSYLTWYQQKQGKSPQLLVYNAKTLAEGVPSRFSGSGSGTQFSLKISSLQPEDFGNYYCQHHYGTRTFGGGTRLEIK")

    print("\n")
    ordered_ids.sort()
    for name_chain in ordered_ids:
        original_seq = sequences[name_chain]
        seq = sequences[name_chain]

        if options.strip_c_term:
            seq = general.strip_right(seq, options.strip_c_term)

        if options.pad_c_term and options.pad_to_n_residues and options.pad_seq_with_stop_codons:
            n_seq = len(seq)
            seq_add = "**"+options.pad_c_term

            total_required = int(options.pad_to_n_residues)

            total_to_pad = total_required - n_seq - 2
            if len(seq_add) < total_to_pad:
                sys.exit("Sequence not long enough!")

            seq_add = seq_add[0:total_to_pad]+"**"
            seq=seq+seq_add

        elif options.pad_c_term and options.pad_to_n_residues:
            n_seq = len(seq)
            seq_add = options.pad_c_term

            total_required = int(options.pad_to_n_residues)

            total_to_pad = total_required - n_seq
            if len(seq_add) < total_to_pad:
                sys.exit("Sequence not long enough!")

            seq_add = seq_add[0:total_to_pad]
            seq=seq+seq_add

        elif options.pad_c_term:
            seq = seq+options.pad_c_term

        if options.three_letter:
            new_seq = []
            for aa in seq:
                new_seq.append(one_to_three(aa))
            seq = " ".join(new_seq)

        if options.format == "basic":
            if options.output_original_seq:
                outlines.append(options.prefix+name_chain+" : "+original_seq)
            outlines.append( options.prefix+name_chain+" : "+seq+"\n")

        elif options.format == "fasta":
            outlines.append("> "+options.prefix+name_chain)
            if options.output_original_seq:
                outlines.append(original_seq+"\n")

            outlines.append(seq+"\n")

        elif igg_type_format:
            name = "IgG_"+name_chain+"_pFUSE"
            if options.output_original_seq:
                outlines.append(str(i)+". "+name+ " "+original_seq+"\n")
            outlines.append(str(i)+". "+name+ " "+seq+"\n")

        elif options.format == "general_order":
            outlines.append(str(i)+". "+"vec"+ " "+seq+"\n")

        i+=1

    #Print always
    for line in outlines:
        print(line)


    #Adding this, even though this whole file needs to be rewritten.  It is more than 10 yrs old at this point...
    if options.format == 'csv':
        df = pandas.DataFrame.from_dict(csv_data)
        outp = "".join(options.outpath.split('.')[0:-1])
        df.to_csv(f'{outp}.csv')

    else:
        OUT = open(options.outpath, "w")
        for line in outlines:
            OUT.write(line + "\n")
        OUT.close()

    #Print any duplicates
    all_sequences = defaultdict(list)
    for id in sequences:
        all_sequences[sequences[id]].append(id)

    for seq in all_sequences:
        if len(all_sequences[seq]) >1:
            print("\nDuplicate sequence in: "+repr(all_sequences[seq]))
            print(">")
            print(seq)

    for id in sequences:
        if sequences[id][1] == 'P':
            print('Proline in second position - run with option (-correct_for_tev) to add G at N-terminus ')
            print(id)

    print("Done.")

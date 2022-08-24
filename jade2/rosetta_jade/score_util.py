import os,sys,re
from collections import defaultdict
from jade2.basic.path import open_file
from jade2.basic.path import get_decoy_name, get_decoy_path
import pandas, json
from jade2.basic.dataframe.util import detect_numeric


def get_dataframe_from_json_csv_pkl(filename: str)-> pandas.DataFrame:
    if filename.split('.')[-1] == "pkl":
        return pandas.read_pickle(filename)
    else:
        return get_dataframe_from_json_or_csv(filename)

def get_dataframe_from_csv(filename: str, test = False) -> pandas.DataFrame:
    print("Reading CSV")
    if test:
        df = pandas.read_csv(filename, nrows = 5000)
    else:
        df = pandas.read_csv(filename)
    df = detect_numeric(df)
    df = create_decoy_path_column(df, filename)
    return df

def get_dataframe_from_json_or_csv(filename: str, test = False)-> pandas.DataFrame:
    if filename.split('.')[-1] == "csv":
        df = get_dataframe_from_csv(filename, test)
    else:
        df = get_dataframe_from_json(filename)
    return df

def get_dataframe_from_json(filename: str) -> pandas.DataFrame:
    """
    Parse a json file and return a dataframe. Can work on any json file.
    If a scorefile, adds decoy path.
    """
    if not os.path.exists(filename):
        sys.exit(filename+" does not exist! Could not parse dataframe!")

    lines = open(filename, 'r').readlines()
    #print(lines)

    #print repr(headerSP)
    decoys = []
    for line in lines:
        line = line.replace("nan", "NaN")
        o = json.loads(line)
            # print o[self.decoy_field_name]
            # print repr(o)
        decoys.append(o)

    df = pandas.DataFrame.from_dict(decoys)
    df = detect_numeric(df)
    #df.to_csv("debugging.csv", sep=",")
    #print(df.columns)

    df["name"]  = os.path.basename(filename)
    df["scorefile"] = os.path.abspath(filename)

    return create_decoy_path_column(df, filename)

def create_decoy_path_column(df: pandas.DataFrame, filename: str) -> pandas.DataFrame:
    def fix_path(original):
        v = get_decoy_path(os.path.dirname(os.path.abspath(filename)) + "/" + original)
        print(v)
        return v

    def maybe_fix_path(original):
        if os.path.exists(original): return original
        else:
            v = os.path.dirname(os.path.abspath(filename)) + "/" + original
            if os.path.exists(v):
                return v
            else:
                v = get_decoy_path(v)
                if os.path.exists(v):
                    return v

        return original

    #Don't overwrite decoy path if it exists.
    if 'decoy_path' in df.columns:
        print("Fixing decoy path")

        df['decoy_path'] = [maybe_fix_path(x) for x in df['decoy_path']]
        if 'decoy' not in df.columns:
            df['decoy'] = [os.path.basename(x) for x in df['decoy_path']]


    elif "decoy" in df.columns:
        if os.path.dirname(filename):

            df["decoy_path"] = [maybe_fix_path(p) for p in df['decoy']]
            print(df["decoy_path"])
            df['decoy_path']= df['decoy_path'].astype(str)
        else:
            df["decoy_path"] = df["decoy"]
            df["decoy_path"] = [get_decoy_path(p) for p in df["decoy_path"]]

    return df

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
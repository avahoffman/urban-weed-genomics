"""
Note! Usage from the command line should be as follows, for example:
  
$> python3 04.5-parse_process_radtags.py --file 'process_radtags.52618524.out'
"""
# import itertools
import argparse
import pandas as pd

def parse_args():
    """
    This function takes the arguments provided at the command line and parses them to use
    below.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file", 
        help="The log output file should go here. E.g.:'process_radtags.52618524.out'", 
        type=str, 
        required=True
    )

    return parser.parse_args()


def generate_stats(filename: str) -> pd.DataFrame:
    """
    This function does everything. It reads the log file from Stacks' process_radtags
    function and then rearranges everything to a nice data frame
    filename: str, refers to the log file that is being looked at.
    returns: DataFrame, containing five columns:
      1. Sublibrary name
      2. Total sequences
      3. Ambiguous barcodes
      4. Low quality read drops
      5. Retained reads
    """

    # Name columns and make an empty dataframe
    cols = ["file", "total", "ambiguousbarcodes", "dropped", "retained"]
    df = pd.DataFrame(columns=cols)

    # Open the file
    Lines = open(filename, "r").readlines()
 
    # Detect lines with filenames and the output results.
    # Strip them and append to the empty dataframe created above
    count = 0
    for line in Lines:
        count += 1
        if "Will attempt to recover" in line:
            filename = Lines[count].strip().strip("Processing file 1 of 1 [").strip("R1.1.fastq.gz]").strip("_")
        if " total sequences" in line:
            total = line
        if " barcode not found drops " in line:
            ambiguousbarcodes = line
        if " low quality read drops " in line:
            dropped = line
        if "retained reads " in line:
            retained = line
            df = df.append(pd.DataFrame([[filename, total, ambiguousbarcodes, dropped, retained]], columns=cols))
            # Remove this print statement after debugging 
            print(df)

    # Split the output results line into actual data
    split_ = df['total'].str.split(' total sequences', expand=True)
    total_ = split_[[0]].rename(columns = {0:'total.sequences'})
    
    split_ = df['ambiguousbarcodes'].str.split(' barcode not found drops ', expand=True)
    ambiguousbarcodes_ = split_[[0]].rename(columns = {0:'ambiguous.barcodes'})
     
    split_ = df['dropped'].str.split(' low quality read drops ', expand=True)
    dropped_ = split_[[0]].rename(columns = {0:'dropped.reads'})
    
    split_ = df['retained'].str.split(' retained reads ', expand=True)
    retained_ = split_[[0]].rename(columns = {0:'retained.reads'})

     # Combine the split columns with the filenames
    out = pd.concat([df['file'], total_, ambiguousbarcodes_, dropped_, retained_], axis=1)
     
    return out
  

def main():
    """
    Wrapper to run the whole thing :)
    """

    args = parse_args()

    print(generate_stats(filename=args.file))


if __name__ == "__main__":
    main()

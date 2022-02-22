"""
Note! Usage from the command line should be as follows, for example:
  
$> python 03.5-parse_clone_filter.py --file 'clone_filter.52618524.out'
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
        help="The log output file should go here. E.g.:'clone_filter.52618524.out'", 
        type=str, 
        required=True
    )

    return parser.parse_args()


def generate_stats(filename: str) -> pd.DataFrame:
    """
    This function does everything. It reads the log file from Stacks' clone_filter
    function and then rearranges everything to a nice data frame
    filename: str, refers to the log file that is being looked at.
    returns: DataFrame, containing three columns:
      1. Sublibrary name
      2. Input number of paired reads
      3. Cleaned up/output number of paired reads.
    """

    # Name columns and make an empty dataframe
    cols = ["file", "stats"]
    df = pd.DataFrame(columns=cols)

    # Open the file
    Lines = open(filename, "r").readlines()
 
    # Detect lines with filenames and the output results.
    # Strip them and append to the empty dataframe created above
    count = 0
    for line in Lines:
      count += 1
      if "Reading data from" in line:
        filename = Lines[count].strip().strip("_R1.fastq.gz and")
      if "pairs of reads input" in line:
        stats = line
        df = df.append(pd.DataFrame([[filename, stats]], columns=cols))

    # Split the output results line into actual data
    split_ = df['stats'].str.split(' pairs of reads input. ', expand=True)
    input_ = split_[[0]].rename(columns = {0:'input'})
    output_ = split_[1].str.split(' pairs of reads output, ', expand=True)[[0]].rename(columns = {0:'output'})

    # Combine the split columns with the filenames
    out = pd.concat([df['file'], input_, output_], axis=1)

    return out
  

def main():
    """
    Wrapper to run the whole thing :)
    """

    args = parse_args()

    print(generate_stats(filename=args.file))


if __name__ == "__main__":
    main()

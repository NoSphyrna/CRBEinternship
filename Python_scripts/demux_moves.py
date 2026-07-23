# This script allows to trim and reverse move tables after the trimmming made by Cutadapt

# TODO Add an 'in place' option with temporary files to avoid multiple folders
# TODO Add more comments and specific comments for the functions
# TODO Parallelise this script to accelerate it (with a ProcessPoolExecutor) see https://realpython.com/python-concurrency/

import argparse
import csv

# We set the buffer of a row to 1 Go to avoid overflow error
_ = csv.field_size_limit(10**9)

import os
import sys

# import pprint
from collections import defaultdict
from glob import glob
from os.path import basename

from Bio.SeqIO.QualityIO import FastqGeneralIterator

# ========================== Parsing of arguments ============================== #
parser = argparse.ArgumentParser(
    description="adapt move tables to adapters trimming by Cutadapt"
)
# Arguments and options of the script
_ = parser.add_argument(
    "input_info_file", type=str, help="Input info file (tsv file)"
)  # the _ = is just to properly discard the result of parser.add_argument (an Action object)
_ = parser.add_argument(
    "input_directory", type=str, help="Input folder with demultiplexed fastq files"
)
_ = parser.add_argument("output_directory", type=str, help="Output directory")


args = parser.parse_args()

info_file = str(args.input_info_file)

# Check the inputs files
if not os.path.isfile(info_file):
    print("\n[python] Error: Info file '", info_file, "' doesn't exist")
    sys.exit(1)
if not info_file.lower().endswith(".tsv"):
    print("\n[python] Error: Input file must have .tsv extension.")
    sys.exit(1)

input_dir = str(args.input_directory)
if not os.path.isdir(input_dir):
    print("\n[python] Error: Fastq dir '", input_dir, "' doesn't exist")
    sys.exit(1)

output_dir = str(args.output_directory)
# If the direcory does not exist, we create it
os.makedirs(output_dir, exist_ok=True)

print("Trimming of mv tables")
# ==================== Parse the info file in a dictionnaire ======================== #

PREFIX_MV = "mv:B:c"
PREFIX_CT = "CT:r:"


def parse_info_file(
    info_file: str, PREFIX_MV: str, PREFIX_CT: str
) -> defaultdict[str, dict[str, int | bool]]:
    # the use of a default dict just handles the case where a key is not found and
    # returns a default value specified by the lambda in that case
    trims = defaultdict(lambda: {"5'": 0, "3'": 0, "rc": False})

    with open(info_file) as file:
        tsv_reader = csv.reader(file, delimiter="\t")  # The info file is a tsv file
        i_id = 0  # The id of the id of the read
        i_check = i_begin = i_end = i_right = i_rc = 0
        row = next(tsv_reader, None)

        while row != None:  # while there are reads to read
            if i_check == 0:  # first row : get the indices
                i_mv, _ = next(
                    (
                        (i, el)
                        for i, el in enumerate(row)
                        if str.startswith(el, PREFIX_MV)
                    ),
                    (0, None),
                )
                if i_mv == 0:
                    print("No mv tables found, skipping trimming of mv tables")
                    sys.exit(0)
                i_check, _ = next(
                    (
                        (i, el)
                        for i, el in enumerate(row)
                        if str.startswith(el, PREFIX_CT)
                    ),
                    (0, None),
                )
                if i_check == 0:
                    print(
                        "Error : No custom field CT:r: found, skipping trimming of mv tables (be sure to add the custom field CT:r: at the end of the header)"
                    )
                    sys.exit(1)
                i_check += 1  # because the check field is the first after the header (so one field after)
                i_begin = i_check + 1  # The id of the begining of the adapter found
                i_end = i_check + 2  # The id of the end of the adapter found
                i_right = (
                    i_check + 5
                )  # The part at the right of the adapter (useful for the 3' adapter)
                i_rc = (
                    i_check + 10
                )  # The boolean wether the read has been reverse complemented

            if int(row[i_check]) != -1:  # if the read contains both adapters

                # Get the infos on the first adapter (5')
                read_id = str(row[i_id])
                # The part trimmed is all bases before the end of the 5' adapter
                n5 = int(row[i_end])
                rc = bool(int(row[i_rc]))

                # The second adpater is on a new row (3')
                row = next(tsv_reader, None)
                if (row == None) or (
                    read_id != str(row[i_id]) or rc != bool(int(row[i_rc]))
                ):
                    print(
                        "Error : Adapters kept but not paired at info fil line : ",
                        tsv_reader.line_num,
                    )
                    sys.exit(1)
                    # might get rid of that later, a bit extreme here but it's for the test
                seq_len = int(row[i_end]) + len(str(row[i_right]))
                # The part trimmed is all bases after the start of the 3' adapter
                n3 = seq_len - int(row[i_begin])
                trims[read_id]["5'"] = n5
                trims[read_id]["3'"] = n3
                trims[read_id]["rc"] = rc

            row = next(tsv_reader, None)
        file.close()
    return trims


trims = parse_info_file(info_file, PREFIX_MV, PREFIX_CT)
# pp = pprint.PrettyPrinter()
print("Info file is parsed")
# pp.pprint(trims)

# ================================ change move tables ================================== #


def trim_mv(
    file: str,
    output_file: str,
    trims: defaultdict[str, dict[str, int | bool]],
    PREFIX_MV: str,
):
    idx_mv = 0
    with open(file) as handle, open(output_file, mode="w") as output_handle:
        for title, sequence, quality in FastqGeneralIterator(handle):
            tags = title.split("\t")
            if idx_mv == 0:
                idx_mv, mv = next(
                    (
                        (i, el)
                        for i, el in enumerate(tags)
                        if str.startswith(el, PREFIX_MV)
                    ),
                    (
                        0,
                        None,
                    ),  # This is the default values returned in case mv is not found
                )
                if mv == None:
                    print("Error : mv field not found in file '", mv, "'")
                    sys.exit(1)
            else:
                mv = tags[idx_mv]
            mv = mv.split(",")
            stride = mv[1]
            mv_table = list(map(int, mv[2:]))
            id = tags[0]
            infos = trims[id]
            n5 = infos["5'"]
            n3 = infos["3'"]
            rc = infos["rc"]
            if rc:
                mv_table.reverse()
            ones = 0
            five_remove = 0
            while ones < n5:
                if mv_table[five_remove]:
                    ones += 1
                five_remove += 1
            ones = 0
            three_remove = 1
            while ones < n3:
                if mv_table[-three_remove]:
                    ones += 1
                three_remove += 1
            three_remove -= 1  # beacause the split of a list is this way :
            # l[a:b] -> keep[a,b[

            mv_table = mv_table[five_remove:-three_remove]
            if len(sequence) != sum(mv_table):
                print(
                    "Error : length of the sequence : ",
                    len(sequence),
                    " isn't equal to length of mv table : ",
                    sum(mv_table),
                    "",
                )
                sys.exit(1)

            tags[idx_mv] = ",".join(
                [PREFIX_MV] + [str(stride)] + list(map(str, mv_table))
            )
            new_title = "\t".join(tags)

            _ = output_handle.write(f"@{new_title}\n{sequence}\n+\n{quality}")
        handle.close()
        output_handle.close()


fastq_files = [f for f in glob(input_dir + "/*.fastq")]
total = len(fastq_files)
print("Nb files to do :", total)
count = 0
for file in fastq_files:
    trim_mv(file, output_dir + basename(file), trims, PREFIX_MV)
    count += 1
    print(count, "/", total, "files treated")

print("Trimming of mv tables is completed")

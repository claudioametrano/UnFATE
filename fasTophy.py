from Bio import AlignIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_fas", default="/dev/stdin")
    parser.add_argument("--out_phy", default="/dev/stdout")

    return parser.parse_args()


def main():

    args = parse_args()

    with open(args.in_fas) as handle:
        records = AlignIO.parse(handle, "fasta")

        with open(args.out_phy, "w") as output_handle:
            AlignIO.write(records, output_handle, "phylip")


if __name__ == "__main__":
    main()

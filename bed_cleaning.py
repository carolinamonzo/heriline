#! /usr/bin/python3.5

# Carolina Monzo; 2018-07-04

def main():

    # Input bed file to clean

    bedfile = str(input("Enter bed file to clean: "))
    with open(bedfile, "r") as fi:
        bed = fi.read().splitlines()

    bedsplit = []

    # Create list of lists to parse the file

    for line in bed:
        bedsplit.append(line.split("\t"))

    sex = ["M", "Y", "X"]
    autosom = []
    for i in range(1, 23):
        autosom.append(str(i))

    chrom = autosom + sex

    good_bed = []

    # Keep only entries where chromosome is autosomal or sexual, not genome patches

    for pos in bedsplit:
        for el in chrom:

            if pos[0] == el:
                good_bed.append("\t".join(pos))


    # Set output name and write clean bed file output


    name = bedfile.split(".")[0]

    with open("./{}_clean.bed".format(name), "a") as f:
        f.write("\n".join(good_bed))


if __name__ == '__main__':
    main()

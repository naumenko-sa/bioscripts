#/usr/bin/env python3

"""
blast.duplex.txt blast.circle.txt threshold_bp
"""

import sys

with open(sys.argv[1], "r") as f:
    l_duplex = f.readlines()

with open(sys.argv[2], "r") as f:
    l_circle = f.readlines()


print("DIST\tDUPLEX_NAME\tDUPLEX_CHR\tDUPLEX_POS\tCIRCLE_NAME\tCIRCLE_CHR\tCIRCLE_POS")

for s_duplex in l_duplex:
    for s_circle in l_circle:
        a_duplex = s_duplex.split("\t")
        a_circle = s_circle.split("\t")

        duplex_name = a_duplex[0]
        duplex_chr = a_duplex[1]
        duplex_pos = a_duplex[6]

        circle_name = a_circle[0]
        circle_chr = a_circle[1]
        circle_pos = a_circle[6]
        
        dist = abs(int(circle_pos) - int(duplex_pos))
        if duplex_chr == circle_chr and dist <= int(sys.argv[3]):
            print(f"{dist}\t{duplex_name}\t{duplex_chr}\t{duplex_pos}\t{circle_name}\t{circle_chr}\t{circle_pos}")


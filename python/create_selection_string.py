"""
Utility script for making residue selection strings for VMD
"""
import os


# residues = range(108, 115)
residues = [35,37,90,38,70,39,71,72,45,14,56,58,91,92,93]

string = "residue "
for i in range(0, len(residues)):
    string = string + str(residues[i]) + " or residue "
string = string[:-12]

os.system("echo '%s' | clip" % string)
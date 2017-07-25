"""
Utility script for making residue selection strings for VMD
"""
import os


# residues = range(108, 115)
residues = [14,35,36,37,38,39,40,45,48,55,56,57,58,69,70,71,72,84,89,90,91,92,93]

string = "residue "
for i in range(0, len(residues)):
    string = string + str(residues[i]) + " or residue "
string = string[:-12]

os.system("echo '%s' | clip" % string)
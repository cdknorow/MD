import os
import sys
sys.path.append("/home/cdknorow/Dropbox/Software/MD/dna_scripts/fitting")
import q6q4_bondorder

for i in sorted(os.listdir('.')):
    if os.path.isdir(i):
        os.chdir(i)
        q6q4_bondorder._run(15.8,5)
        os.chdir('../')

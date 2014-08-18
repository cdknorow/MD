import os
import sys
sys.path.append("/home/cdknorow/Dropbox/Software/MD/dna_scripts/fitting")
import q6q4_bondorder

for i in sorted(os.listdir('.')):
    if os.path.isdir(i):
        os.chdir(i)
        if not os.path.isfile('q6q4_all_rcut7.00.xyz'):
            try:
                #sp4
                #q6q4_bondorder._run(11.75,10)
                #sp8
                #q6q4_bondorder._run(15.80,5)
                #sp12
                #q6q4_bondorder._run(17.2,5)
                #sp2
                #q6q4_bondorder._run(8.5,5)
                #sp0
                #q6q4_bondorder._run(5.85,5)
                ## SOFT ##
                #sp4
                #q6q4_bondorder._run(11.5,10)
                #sp8
                q6q4_bondorder._run(14,5)
                #sp12
                #q6q4_bondorder._run(,5)
                #sp2
                #q6q4_bondorder._run(,5)
                #sp0
                #q6q4_bondorder._run(5.85,5)
            except:
                print 'Skipping',i
        else:
            print 'q6q4 already exists'
        os.chdir('../')

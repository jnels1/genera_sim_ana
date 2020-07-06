import pickle
import sys

sys.path.append('/afs/slac.stanford.edu/u/dm/jnels1/general_sim_ana')
import general_tools

def main():

    infile = sys.argv[1]
    outfile = sys.argv[2]

    n_det = 24

    skimmed_data = general_tools.skim_all(infile, n_det)

    print('writing hits from '+str(infile)+' to '+str(outfile))
    with open(outfile, 'wb') as f:
        pickle.dump( skimmed_data, f )

    return;

if __name__=="__main__":
    main()

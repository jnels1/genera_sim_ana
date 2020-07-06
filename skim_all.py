import pickle
import sys
import os
import datetime

sys.path.append('/afs/slac.stanford.edu/u/dm/jnels1/general_sim_ana')
import general_tools

def main():
    # Get date for labeling the savefile
    x = datetime.datetime.now()
    date_str = x.strftime("%d")+x.strftime("%m")+x.strftime("%y")
    
    #indir = sys.argv[1]
    #outdir = sys.argv[2]
    #save_name = sys.argv[3]
    indir = '/nfs/slac/g/cdms/u05/env_sims/replay_data/biased/no_shield/data/'
    outdir = '/nfs/slac/g/cdms/u05/env_sims/replay_data/biased/no_shield/combined/'
    save_name = 'fixed_rethrow_combo' 

    n_det = 24

    combo_dict=dict()

    file_list = os.listdir(indir)
    for name in file_list:

        skimmed_data = general_tools.skim_all(indir+name, n_det)
        print('writing hits from '+str(name))
        combo_dict[name]=skimmed_data


    with open(outdir+save_name+'_'+date_str+'.p', 'wb') as f:
        pickle.dump( combo_dict, f )

    return;

if __name__=="__main__":
    main()

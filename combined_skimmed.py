import pickle
import sys
import os

sys.path.append('/afs/slac.stanford.edu/u/dm/jnels1/general_sim_ana')
import general_tools

def combine_util( skimmed_data_dir, outname  ):
    
    file_list = os.listdir(skimmed_data_dir)
    save_dir = '/nfs/slac/g/cdms/u05/env_sims/cavern_gammas/c_stem_combined/'
    combo_dict=dict()
    
    for name in file_list:
        with open(skimmed_data_dir+name, 'rb') as f:
            temp_data = pickle.load(f)
        combo_dict[name]=temp_data

    with open(save_dir+outname+'.p', 'wb') as f:
        pickle.dump(combo_dict, f)

    return;

def main():
    
    basedir = '/nfs/slac/g/cdms/u05/env_sims/cavern_gammas/'
    # Combine the background data
    combine_util( basedir+'c_stem_skimmed/', 'bg_data_combo' )
    
    # Combine the replay data
    combine_util( basedir+'c_replay_skimmed/', 'replay_data_combo' )

    return;

if __name__=="__main__":
    main()

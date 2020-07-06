import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import uproot
import pickle
import time
import sys
import os

def get_cylinder( data_frame, vol_name ):
    
    tol = 0.003 # tolerence for the edge. This choice is 0.5 the thickness of the can bottom
    r3 = np.sqrt( np.power(data_frame['X3'], 2) + np.power(data_frame['Y3'], 2) )

    theta3 = np.arctan2( data_frame['Y3'], data_frame['X3'] )

    volume_dims = { 'mc':{'r0':0.20955, 'z_top':-0.02370, 'z_bot':-0.464934},
                    'poly':{'r0':1.36246, 'z_top':1.1298133, 'z_bot':-1.13436} }

    is_in = data_frame['InOut'] == 1
    is_out = data_frame['InOut'] == 2
    is_through = data_frame['InOut'] == 3

    vol_map = {'poly':b'Scorers/Shield/InnerPoly', 'mc':b'Scorers/Vessel/MixingChamber'}
    vol_oi = data_frame['VolName']==vol_map[vol_name]

    #print(np.count_nonzero(is_in), np.count_nonzero(is_out), np.count_nonzero(is_through) )


    z_top = volume_dims[vol_name]['z_top']
    z_bot = volume_dims[vol_name]['z_bot']
    r0 = volume_dims[vol_name]['r0']

    inner_region = ((r3<=r0-tol))&((data_frame['Z3']<=z_top-tol))&((data_frame['Z3']>=z_bot+tol))
    outer_region = ~inner_region

    c_side = data_frame['X3'] < 0
    e_side = data_frame['X3'] > 0
  
    cyl_gammas = data_frame.loc[(is_out|is_through)&inner_region&vol_oi]
    cyl_gammas['r'] = r3[(is_out|is_through)&inner_region&vol_oi]
    cyl_gammas['theta'] = theta3[(is_out|is_through)&inner_region&vol_oi]
    
    return cyl_gammas;

def get_shield_hits(root_file):
    
    dim_dict = pd.read_csv( 'cylinder_dim_ref.csv' )

    # load the scorer data from the root file
    score_obj = uproot.open(root_file)['G4SimDir/mcScorer']
    score_data = score_obj.pandas.df("*")
    score_data['File'] = root_file

    event_obj = uproot.open(root_file)['G4SimDir/mcevent']
    event_data = event_obj.pandas.df("*")
    prims = len(event_data['EventNum'])

    vol_map = {'mc':b'Scorers/Vessel/MixingChamber',
               'poly':b'Scorers/Shield/InnerPoly'}
    
    shield_dict = {}
    
    for i, name in enumerate(list(vol_map)):
        dim_index = np.where(dim_dict['volume']==name)[0]
        shield_dict[name] = get_cylinder( score_data, name )
    
    shield_dict['primaries']=prims
        
    return shield_dict;
# Finished Functions

def grab_zip_data(root_file, zip_num):
    zip_obj = uproot.open(root_file)['G4SimDir/mczip'+str(zip_num)]
    zip_data = zip_obj.pandas.df("*")
    zip_data['File'] = root_file
    not_empty = zip_data['Empty']==0
    return zip_data.loc[not_empty];

def grab_prim_data(root_file):
    prim_obj = uproot.open(root_file)['G4SimDir/mcprimary']
    prim_data = prim_obj.pandas.df("*")
    return prim_data;

def grab_decay_data(root_file):
    dec_obj = uproot.open(root_file)['G4SimDir/mcDecays']
    dec_data = dec_obj.pandas.df("*")
    return dec_data;

def group_hits( zip_data, time_window):  
    # first, identify time gaps above the defined window
    time_break = np.diff(zip_data['Time1'])>time_window
    break_index = np.where(time_break)[0]

    # assign a 'group index' which increases by 1 each break
    time_group = np.zeros(len(zip_data))
    for i in range(np.count_nonzero(time_break)):
        # There's probably a cleaner way to do this
        temp_group = np.ones(len(zip_data))
        temp_group[:break_index[i]+1] = 0
        time_group = time_group + temp_group

        zip_data['TimeGroup'] = time_group
        
        # Now, use the pd.DataFrame.groupby() method to sort the data by event number and time group
        # Each of these should represent an analysis event
        
        grouped_data = zip_data.groupby(['TimeGroup', 'EventNum'])
    
    return grouped_data;

def make_rq( analysis_event ):
    # take an analysis event and extract some "RQs"
    #print(list(analysis_event))
    time_start = np.amin( analysis_event['Time1'] )
    # classify the amount of energy deposited by electron recoils and nuclear recoils
    is_n = (( analysis_event['PType'] >= 1e4) | (analysis_event['PType'] == 2112))
    is_e = (( analysis_event['PType'] < 1e4 ) & (analysis_event['PType'] != 2112)) # should just be ~is_nr, but i'll write it out for clarity

    nr_tot = np.sum( analysis_event.loc[ is_n, 'Edep'] )
    er_tot = np.sum( analysis_event.loc[ is_e, 'Edep'] )    
    
    rq = { 'Time':time_start,
           'ETotal':er_tot+nr_tot,
           'ER_part':er_tot/(er_tot+nr_tot),
           'NR_part':nr_tot/(er_tot+nr_tot),
           'DetNum':analysis_event['DetNum'].iloc[0],
           'EventNum':analysis_event['EventNum'].iloc[0],
           'File':analysis_event['File'].iloc[0]
           }
                                                                                                                                                                                                                                        
    # check if there's decay information
    # Right now, I only have info for the U238 chain implemented 
    decay_tree = { '238092':'U238',
                   '234090':'Th234',
                   '234091':'Pa234',
                   '234092':'U234',
                   '230090':'Th230',
                   '226088':'Ra226',
                   '222086':'Rn222',
                   '218084':'Po218',
                   '218085':'At218',
                   '214086':'Rn218',
                   '214082':'Pb214',
                   '214083':'Bi214',
                   '214084':'Po214',
                   '210081':'Tl210',
                   '210082':'Pb210',
                   '210083':'Bi210',
                   '210084':'Po210',
                   '206080':'Hg206',
                   '206081':'Tl206',
                   '206082':'Pb206',
                   '0':'NA' }
    
    if 'decayAncestor.PType' in list(analysis_event):
        unq_events = np.unique(analysis_event['decayAncestor.PType'])
        if len(unq_events)==1:
            # We only see one decay parent in the hit, which makes our job easy
            first_parent = int(analysis_event['decayAncestor.PType'].iloc[0])
            is_mixed=False
        elif len(unq_events)==2:
            # Two decay parents, infer whether they're related by alpha or beta decay
            A = np.floor(unq_events/100000.)
            Z = np.mod(unq_events, 100000.)
            if A[0]==A[1]:
                # beta decay
                first_parent = int(unq_events[np.argmin(Z)])
            else:
                # assume it's an alpha decay
                first_parent = int(unq_events[np.argmax(unq_events)])
                is_mixed=True
        else:
            # anything more... don't bother
            first_parent = 0
            is_mixed=True
            
        rq['ParentList']=np.unique(analysis_event['decayAncestor.PType'])
        rq['Decay']=decay_tree[str(first_parent)]
        rq['MixedDecay']=is_mixed
    
    return rq;

def get_one_zip_rq( root_file, zip_num ):
    hit_list = []
    
    zip_data = grab_zip_data(root_file, zip_num)
    
    if len(zip_data)>0:
        hits_obj = group_hits(zip_data, 10000)
        print('processing '+str(len(hits_obj))+' in det. '+str(zip_num))
        for name, hit in hits_obj:
            p_count = make_rq(hit)
            hit_list.append( p_count )
        
    return pd.DataFrame(hit_list);

def get_one_file_rq( root_file, n_zips ):
    
    # Inputs:
    #    root_file: just the path to the file we want to analyze
    #    n_zips: the number of detectors in the simulation
    
    zip_nums = np.arange(1, n_zips+1)
    
    # now we can loop over all the detectors and collect the hits
    
    all_hits=dict()
    print('Found '+str(zip_nums)+' detectors')
    for i in zip_nums:
        #print('getting data from det. '+str(i))
        all_hits['Det'+str(i)] = get_one_zip_rq(root_file, i)
    return all_hits;

def check_mult( rq_data, mult_window ):
    
    # for each detector and event, check if there are any other hits in another detector within
    # the specified time window
    
    # Note that techinically this only indicates a potential multiple. For it to be a true multiple, both hits need
    # to be above the detector threshold. I don't want to specify a threshold now for generality. So just keep in
    # mind that this is only half the boolean for a real analysis multiple
    
    rq_mult = rq_data.copy()
    
    # Loop over each detector
    for det_name in rq_data:
        # most events are singles. Initialize to 'true'
        rq_mult[det_name]['Single']=True
        # Need to loop over each event in the detector
        for i, event in enumerate( rq_data[det_name]['EventNum'] ):
            # And loop over each other detector
            for co_det in rq_data:
                # Cant have multiples in the same detector (does this make sense?)
                if det_name==co_det:
                    continue
                # check if this detector has any hits in this event
                this_event = rq_data[co_det]['EventNum']==event
                # if there aren't any hits, move on
                if np.count_nonzero(this_event)==0:
                    continue
                else:
                    # Get the time difference between the hit of interest and all hits in the same event in the
                    # corresponding detector
                    time_delta = np.abs(rq_data[det_name]['Time'].iloc[i] - rq_data[co_det]['Time'].loc[this_event])
                    # if any of this is less than the threshold, thats a multiple, baby!
                    if np.count_nonzero(time_delta<mult_window):
                        rq_mult[det_name]['Single'].iloc[i]=False
                        
    return rq_mult;

def skim_all(root_file, n_det):

    shield_hits = get_shield_hits(root_file)
    try:
        det_hits = get_one_file_rq(root_file, n_det)
    except:
        det_hits = pd.DataFrame()
    ret_dict = {'shield':shield_hits, 'zips':det_hits}

    return ret_dict;

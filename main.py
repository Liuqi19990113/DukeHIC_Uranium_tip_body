import numpy as np
import h5py
import obs
import sys
import os

eta_std = [-1, 1]
eta_sub1 = [-1,-0.3]
eta_sub2 = [0.3, 1]
y_cut = [-0.1, 0.1]
pt_cut = [0.2, 2]

def Go(hdf5_file_list,mult_cut):
    results_list = []
    k = 0
    j = 0
    m = 0
    for file in hdf5_file_list:
        k+=1
        print('Begin {}th file'.format(k))
        with h5py.File(file, 'r') as f:
            for ev_name in f.keys():
                nsample = f[ev_name]['sample'][-1]
                sample = f[ev_name]['sample']
                charge = f[ev_name]['charge']
                phi = f[ev_name]['phi']
                eta = f[ev_name]['eta']
                pt = f[ev_name]['pT']
                y = f[ev_name]['y']
                spin_a, spin_b, tilt_a, tilt_b = f[ev_name].attrs['spin_a'],f[ev_name].attrs['spin_b'],\
                                                 f[ev_name].attrs['tilt_a'],f[ev_name].attrs['tilt_b']
                for i in range(1,nsample+1):
                    m+=1
#                    print('Begin {}th event'.format(j))
                    mult_ref = (phi[(sample == i) & (charge != 0) & (eta > -0.5) & 
                            (eta < 0.5) ] ).size  #You can change the condition.
                    if mult_ref < mult_cut:
                        continue
                    else:
                        j+=1
                        print('{} events added'.format(j))
                        phi_std_ref = phi[(sample == i) & (charge != 0) & (eta > eta_std[0]) & 
                            (eta < eta_std[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                        phi_sub1 = phi[(sample == i) & (charge != 0) & (eta > eta_sub1[0]) & 
                            (eta < eta_sub1[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                        phi_sub2 = phi[(sample == i) & (charge != 0) & (eta > eta_sub2[0]) & 
                            (eta < eta_sub2[1]) & (pt > pt_cut[0]) & (pt < pt_cut[1])]  #You can change the condition.
                        pT_array = pt[(sample == i) & (charge != 0) & (y > y_cut[0]) & (y < y_cut[1])]
                        ave2_2,ave4_2,ave2_3,ave4_3,ave2_4,ave4_4,meanpT = obs.runobs(phi_std_ref,phi_sub1,phi_sub2,pT_array)
                        results_list.append(np.array([mult_ref,ave2_2,ave4_2,ave2_3,ave4_3,ave2_4,ave4_4,meanpT,spin_a, spin_b, tilt_a, tilt_b]))
    print('total: {}, add: {}'.format(m,j))
    return np.array(results_list)


if __name__ == '__main__':
    dicname = sys.argv[1:]
    final_list = []
    for dic in dicname:
        tem=os.listdir(dic)
        new_tem=[os.path.join(dic,name) for name in tem]
        final_list+=new_tem
    results=Go(final_list,695)
    flc_pT = (results[:,7]-np.mean(results[:,7]))/np.mean(results[:,7])
    results=np.insert(results,8,flc_pT,axis=1)
    np.savez('results.npz',results)
    
        
    



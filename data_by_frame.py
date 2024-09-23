import mdtraj as md
import numpy as np
import pandas as pd
from multiprocessing import Pool
import scipy
import math
import argparse
import os
import time

start = time.time()

SIDE_CHAINS_LIST = ["GLY", "ALA", "VAL", "LEU", "ILE", "ASP", "GLU", "ASN", "GLN", "PRO",
                    "PHE", "TRP", "LYS", "CYS", "MET", "TRY", "ARG", "HIS", "SER", "THR"]
ELEMENT_LIST = ["oxygen", "nitrogen", "sulfur"]
DISTANCE_CUTOFF = 3.5
ANGLE_CUTOFF_RAD = 0.523599
PATH = os.getcwd()

# file input
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Generate interaction information file.")
    requirements = parser.add_argument_group("required arguments")
    requirements.add_argument("-top", dest="topology", help="input file, topology file", type=str, required=True)
    requirements.add_argument("-traj", dest="trajectory", help="input file, trajectory file", type=str, required=True)
    requirements.add_argument("-cc", dest="clustercenter", help="input file, clustercenter file", type=str, required=True)
    requirements.add_argument("-summary", dest="summary", help="input file, hsa_hsa_summary (txt or csv)", type=str, required=True)
    requirements.add_argument("-sitewat", dest="site_water", help="input file, site water csv (optional)", type=str, required=False)
    requirements.add_argument("-hs", dest="hs_site", help="hydration site number to make the movie", type=int, required=True)
    requirements.add_argument("-nframes", dest="num_frames", help="the number of frames to be loaded", type=int, required=True)
    requirements.add_argument("-pref", dest="prefix", help="output file prefix", type=str,required=True)
    args = parser.parse_args()

    top_file = args.topology
    trj_file = args.trajectory
    cc = args.clustercenter

    sum = args.summary
    site_wat = args.site_water
    hs = args.hs_site
    nframes = args.num_frames
    pref = args.prefix

top = md.load_topology(top_file)
# trj = md.load(trj_file, top=top)
trj_water = top.select("water and name O")
first_water = int(trj_water[0])
center = md.load(cc)
num_hs_water = center.n_atoms
len_of_traj = range(nframes)
center_crds = center.xyz[0, :, :]*10

########## Part 1: site water

def site_water(i):
    t = md.load_frame(trj_file, i, top=top)
    wat_crd = t.xyz[0, trj_water, :] * 10
    distance = scipy.spatial.distance.cdist(wat_crd, center_crds, "euclidean")
    d = np.where(distance <= 1.0)
    df_site = pd.DataFrame(columns=range(0, num_hs_water), index=range(1))
    for k in range(len(d[0])):
        a = d[0][k]
        b = d[1][k]
        df_site.iloc[0][b] = int(a * 4 + first_water)

    return df_site

if site_wat:
    site_water_data = pd.read_csv(site_wat, index_col="Unnamed: 0")
else:
    p = Pool(6)
    result = p.map(site_water, list(len_of_traj))
    site_water_data = pd.concat(result)
    site_water_data.index = len_of_traj
    site_water_data.to_csv(pref + "_site_water.csv", ",")
    p.close()
    p.join()

t_mid = time.time()
mid_time = round(t_mid - start, 2)
print("Hydration site water processing time: ", mid_time, " sec")

########## Part 2: read hsa summary, side chains making h-bond
#### If the summary file format is csv, run read_csv and process_df.
#### If the summary file format is txt, run read_txt and process_df.

def read_csv(hs_summary_file):
    df_csv = pd.read_csv(hs_summary_file, index_col="index")

    return df_csv

def process_df(df_csv):
    df = df_csv
    sol_acc = df.iloc[:, -2]
    sol_don = df.iloc[:, -1]
    sol_acc_list = []
    sol_don_list = []

    for sola in sol_acc:
        each_alist = []
        if sola != "NONE":
            # res_num_list = []
            acc_list = sola.split(',')
            for acc in acc_list:
                res_a = int(acc.split('-')[0][3:])
                atom_a = acc.split('-')[1]
                atom_a_idx = top.residue(res_a).atom(atom_a).index
                each_alist.append([str(top.residue(res_a)), atom_a_idx])
        sol_acc_list.append(each_alist)

    for sold in sol_don:
        each_dlist = []
        if sold != "NONE":
            res_num_list = []
            don_list = sold.split(',')
            for don in don_list:
                res_d = int(don.split('-')[0][3:])
                atom_d = don.split('-')[1]
                atom_d_idx = top.residue(res_d).atom(atom_d).index
                each_dlist.append([str(top.residue(res_d)), atom_d_idx])
        sol_don_list.append(each_dlist)

    return sol_acc_list, sol_don_list


def read_txt(hs_summary_file):
    if hs_summary_file[-3:] != "csv":
        with open(hs_summary_file, "r") as file_txt:
            lines_txt = file_txt.readlines()
            num_cluster = len(lines_txt) - 1
            index_txt = lines_txt[0]
            to_columns = index_txt.split()
            df_txt = pd.DataFrame(columns=to_columns)
            for nth in range(num_cluster):
                each_line = lines_txt[nth+1]
                data_txt = each_line.split()
                if len(data_txt) == 28:
                    acc_sw = float(data_txt[-3])
                    don_sw = float(data_txt[-2])
                    if acc_sw == 0.0:
                        data_txt.append("NONE")
                    elif don_sw == 0.0:
                        sol_don = data_txt[-1]
                        data_txt = data_txt[:-1]
                        data_txt.append("NONE")
                        data_txt.append(sol_don)
                elif len(data_txt) == 27:
                    num_to_nan = 29-len(data_txt)
                    for nan in range(num_to_nan):
                        data_txt.append("NONE")
                df_txt.loc[str(nth)] = data_txt

        df_txt.to_csv(hs_summary_file[:-3] + "csv", ",", index=False)

    return df_txt

########## Part 3: find frames making h-bond
d=0
def check_criteria(frame, don_at, acc_at):
    global d
    # global a
    # dd = top.atom(don_at).residue
    # dda = top.atom(don_at).name
    # ddah1 = top.atom(don_at+1).name
    # ddah2 = top.atom(don_at+2).name
    # aa = top.atom(acc_at).residue
    # aaa = top.atom(acc_at).name
    # print(don_at, acc_at)
    # print(dd, dda, ddah1, ddah2, aa, aaa)
    don_crd = [frame.xyz[0, don_at, :] * 10]
    acc_crd = [frame.xyz[0, acc_at, :] * 10]
    dist = scipy.spatial.distance.cdist(don_crd, acc_crd, "euclidean")
    if dist <= 3.6:
        angle_triplet_H1 = [[acc_at, don_at, don_at+1]]
        angle_H1 = md.compute_angles(frame, angle_triplet_H1)
        angle_triplet_H2 = [[acc_at, don_at, don_at+2]]
        angle_H2 = md.compute_angles(frame, angle_triplet_H2)
        if angle_H1 <= ANGLE_CUTOFF_RAD:
            d+=1
            return don_at
        elif angle_H2 <= ANGLE_CUTOFF_RAD:
            d+=1
            return don_at
        else:
            return "NoHB"
    else:
        return "NoHB"


def frames_hbond(frame_i, site_waters, solute_acc, solute_don):
    df_frame = pd.DataFrame(columns=["frame", "site_wat_resi", "site_wat_idx", "don_resi", "don_index", "acc_resi", "acc_index"], index=range(1))
    df_frame.iloc[0][0] = frame_i                       # put every frame
    frame = md.load_frame(trj_file, frame_i, top=top)

    if math.isnan(site_waters[frame_i]) is False:       # if there is a site water
        wat_O = int(site_waters[frame_i])
        wat_res = top.atom(wat_O).residue.index
        df_frame.iloc[0][1] = wat_res               # put site water resi if there is a site water
        df_frame.iloc[0][2] = wat_O                 # put water oxygen idx
        if len(solute_acc) != 0:
            acc_res_list = []
            acc_at_n_list = []
            for acc_pair in solute_acc:
                acc_res = acc_pair[0]
                acc_at_n = acc_pair[1]
                hbond = check_criteria(frame, wat_O, acc_at_n)
                # print(hbond)
                if hbond != "NoHB":                     # if water donate, solute accepts
                    acc_res_list.append(acc_res)
                    acc_at_n_list.append(acc_at_n)

            # df_frame.iloc[0][2] = wat_O + 1     # put water oxygen idx
            df_frame.iloc[0][5] = ', '.join(acc_res_list)   # put solute acceptor resi
            df_frame.iloc[0][6] = ', '.join([str(a) for a in acc_at_n_list])   # put solute acceptor idx


        if len(solute_don) != 0:
            don_res_list = []
            don_at_n_list = []
            for don_pair in solute_don:
                don_res = don_pair[0]
                don_at_n = don_pair[1]
                hbond = check_criteria(frame, don_at_n, wat_O)
                if hbond != "NoHB":                     # if water accepts, solute donates
                    don_res_list.append(don_res)
                    don_at_n_list.append(don_at_n)
            df_frame.iloc[0][3] = ', '.join(don_res_list)    # put solute donor resi
            df_frame.iloc[0][4] = ', '.join([str(a) for a in don_at_n_list])   # put solute donor idx


    return df_frame


########## Part 4

if sum[-3:] == "txt":
    txt_read = read_txt(sum)
    acc_data, don_data = process_df(txt_read)
elif sum[-3:] == "csv":
    csv_read = read_csv(sum)
    acc_data, don_data = process_df(csv_read)

def run_fcn(i):
    site_wat = site_water_data.iloc[:, hs]
    result_df = frames_hbond(i, site_wat, acc_data[hs], don_data[hs])
    return result_df

p = Pool(6)
result = p.map(run_fcn, list(len_of_traj))
output = pd.concat(result)
output.index = len_of_traj
output.to_csv(pref + "_hs" + str(hs) + "_hbond_info.csv", ",")
p.close()
p.join()

end = time.time()
total_time = round(end - start, 3)
print("Total time passed: ", total_time, " sec")
print(d)
# to csv
# analyze
# generate pml file
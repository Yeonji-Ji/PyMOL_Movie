import pandas as pd
import numpy as np
import argparse
import math
import pymol

cartoon_color = "grey"
cartoon_transparency = "0.8"
sphere_scale = "1.0"
sphere_transparency = "0.8"

def settings(top, cc, hs):
    pymol.cmd.bg_color("white")
    pymol.cmd.set("depth_cue", 0)
    pymol.cmd.load(top, "top")
    pymol.cmd.load(cc, "cent")
    pymol.cmd.hide("everything", "cent")
    pymol.cmd.select("site", "cent and resi " + str(hs))
    # cmd.show("spheres", "cent and site")
    # cmd.set("sphere_transparency", sphere_transparency, "cent and site")
    # cmd.set("sphere_scale", sphere_scale, "cent and site")
    return

def traj(traj, start, end):
    pymol.cmd.load_traj(traj, "top", start=start+1, stop=end+1)
    pymol.cmd.hide("everything", "top and solvent")
    pymol.cmd.set("cartoon_color", cartoon_color, "top")
    pymol.cmd.set("cartoon_transparency", cartoon_transparency, "top")
    return

def scenemaker(input_file, start, end):
    df = pd.read_csv(input_file, index_col="frame")
    df_site_wat_res = df.iloc[start:end+1, 1]
    df_site_wat_idx = df.iloc[start:end+1, 2]
    df_don_resi = df.iloc[start:end+1, 3]
    df_don_resi_at = df.iloc[start:end+1, 4]
    df_acc_resi = df.iloc[start:end+1, 5]
    df_acc_resi_at = df.iloc[start:end+1, 6]


    for f in range(1, end-start+2):
        if math.isnan(df_site_wat_res.iloc[f-1]) is False:
        # if type(df_site_wat_res.iloc[f-1]) == str:
            wat_resi = "top and resi " + str(int(df_site_wat_res.iloc[f-1])+1)
            # print(wat_resi)
            # if math.isnan(df_don_resi.iloc[f-1]) is False:
            if type(df_don_resi.iloc[f-1]) == str:
                don_resi_obj = '+'.join([str(int(r[3:])+1) for r in df_don_resi.iloc[f-1].split(", ")]) ## (e.g. GLU154) Extract "154" and +1 for pymol and combine with "+" if multiple
                sol_don_resi = "top and resi " + don_resi_obj
                # print(don_resi_obj)
                if type(df_don_resi_at.iloc[f-1]) == str:
                    don_resi_at_obj = '+'.join([str(int(a)+1) for a in df_don_resi_at.iloc[f-1].split(", ")])
                else:
                    don_resi_at_obj = str(int(df_don_resi_at.iloc[f - 1]) + 1)
                sol_don_idx = "top and index " + don_resi_at_obj
            # elif (type(df_don_resi_at.iloc[f-1]) == float) and (math.isnan(df_don_resi_at.iloc[f-1]) is False):
            #     don_resi_at_obj = str(int(df_don_resi_at.iloc[f-1]) + 1)
            #     sol_don_idx = "top and index " + don_resi_at_obj
            # if math.isnan(df_acc_resi.iloc[f-1]) is False:
            if type(df_acc_resi.iloc[f-1]) == str:
                acc_resi_obj = '+'.join([str(int(r[3:])+1) for r in df_acc_resi.iloc[f-1].split(", ")])
                sol_acc_resi = "top and resi " + acc_resi_obj
                if type(df_acc_resi_at.iloc[f - 1]) == str:
            # if type(df_acc_resi_at.iloc[f-1]) == str or float:
                    acc_resi_at_obj = '+'.join([str(int(a)+1) for a in df_acc_resi_at.iloc[f-1].split(", ")])
                else:
                    acc_resi_at_obj = str(int(df_acc_resi_at.iloc[f - 1]) + 1)
                sol_acc_idx = "top and index " + acc_resi_at_obj
            #     sol_acc_idx = "top and index " + acc_resi_at_obj
            # elif (type(df_acc_resi_at.iloc[f-1]) == float) and (math.isnan(df_acc_resi_at.iloc[f-1]) is False):
            #     acc_resi_at_obj = str(int(df_acc_resi_at.iloc[f - 1]) + 1)
            #     sol_acc_idx = "top and index " + acc_resi_at_obj

            pymol.cmd.frame(f)
            pymol.preset.ball_and_stick(selection=wat_resi + " and not name EPW", mode=1)
            pymol.cmd.unbond(wat_resi + " and name H1", wat_resi + " and name H2")
            # if math.isnan(df_don_resi.iloc[f-1]) is False:
            if type(df_don_resi.iloc[f - 1]) == str:
                pymol.cmd.show("sticks", sol_don_resi)
                pymol.cmd.distance("d_acc" + str(f), wat_resi + " and element O", sol_don_idx)
                pymol.cmd.hide("label", "d_acc" + str(f))
                pymol.cmd.color("tv_red", "d_acc" + str(f))
            # if math.isnan(df_acc_resi.iloc[f-1]) is False:
            if type(df_acc_resi.iloc[f - 1]) == str:
                pymol.cmd.show("sticks", sol_acc_resi)
                pymol.cmd.distance("d_don" + str(f), wat_resi + " and element O", sol_acc_idx)
                pymol.cmd.hide("label", "d_don" + str(f))
                pymol.cmd.color("tv_blue", "d_don" + str(f))
            pymol.cmd.orient("site")
            pymol.cmd.scene("scene_" + str(f), "store")
            pymol.cmd.hide("everything", wat_resi)
            pymol.cmd.disable("d_*" + str(f))
        else:
            pymol.cmd.scene("scene_" + str(f), "store")



    return

def scene_view(num_frames, view, resi=None):
    for f in range(1, num_frames+1):
        pymol.cmd.scene("scene_"+str(f))
        pymol.cmd.set_view(view)
        pymol.cmd.show("spheres", "cent and site")
        pymol.cmd.set("sphere_scale", sphere_scale, "cent and site")
        pymol.cmd.set("sphere_transparency", sphere_transparency, "cent and site")
        if resi is not None:
            pymol.cmd.show("sticks", "top and resi "+str(resi))
        pymol.cmd.scene("scene_" + str(f), "update")

    return

def make_movie(num_frames):
    pymol.cmd.mset('1 x7')
    pymol.cmd.scene("scene_1", animate=0)
    pymol.cmd.mview("store", 1, 7)
    for f in range(2, num_frames+1):
        pymol.cmd.madd(str(f) + ' x7')
        pymol.cmd.scene("scene_" + str(f), animate=0)
        pymol.cmd.mview("store", (f-1)*7+1, f*7)

    pymol.cmd.mview("interpolate")
    pymol.cmd.mview("smooth")
    pymol.cmd.set("movie_loop", 0)
    return

# load_files = "load " + top + ", top\n load " + cc + ", cent\n load_traj " + traj + ", top, start=" + start + ", end=" + end + "\n"
# represent_vis = "hide everything, top and solvent\n hide everything, cent\n show spheres, " + center + "\n"
# represent_detail = "set cartoon_color, " + cartoon_color + ", top\n set cartoon_transparency, " + cartoon_transparency + ", top\n set sphere_scale, " + sphere_scale + ", " + center + "\n"
# bg_setting = "set bg_rgb=[1, 1, 1]\n set depth_cue, 0\n"

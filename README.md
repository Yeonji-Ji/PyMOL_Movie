1.	“Data_by_frame.py” script provides the protein – water interaction of a hydration site by each frame. This step is to find the frames worth making movie.

e.g.) python   data_by_frame.py   -top   top.prmtop   -traj   traj.nc  -cc clustercenterfile.pdb    -summary    hsa_hsa_summary.csv    -hs  5    -pref  pref

Arguments \n
    top) top.ptmtop\n
    traj) trajectory.nc\n
    cc) clustercenterfile.pdb\n
    summary) hsa_hsa_summary.txt (or csv if you have)\n
    hs) hydration site number to make a movie\n
    pref) prefix for the data file\n
    
Optional \n
    sitewat) hydration_site_water_file.csv (if you have. If you don’t have it, the script automatically generates the csv file, so you can use the file for the next run)

2.	Open PyMOL.
  
3.	In the command line, write “run write_pymol_movie.py”
** “write_pymol_movie.py” has five functions; settings, traj, scenemaker, scene_view, make_movie.
  	
4.	Run settings(top.prmtop, clustercenterfile.pdb, hydration site number to see) in the command line to set up the files and visualization state.

e.g.) settings(“akt1.prmtop”, “clustercenterfile.pdb”, 5)

5.	Run traj(traj.nc, starting frame, end frame) to load the frames of trajectory required for the movie.

e.g.) traj(“traj.nc”, 25, 59)

6.	Run scenemaker(output file from “data_by_frame.py” run, starting frame, end frame) to generate scenes of the movie. Each scene is one frame loaded and shows you the h-bonding interaction of the hydration site water with the side chain.

e.g.) scenemaker(“output.csv”, 25, 59)

8.	Orient the molecule for the right view by running “my_view=cmd.get_view()”.

9.	Run scene_view(number of frames (or scenes), view, resi) to make the view of scenes the same. And if there is specific residues to show in the whole movie, put them in “resi”.

e.g.) scene_view (35, my_view, 92+84)

10.	Run make_movie(number of frames (or scene)) to generate a movie.

e.g.) make_movie(35)

11.	Save the movie. 

File – Export Movie as.. - MPEG...


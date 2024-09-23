import pymol


def make_movie(num_frames):
    pymol.cmd.mset('1 x7')
    pymol.cmd.scene("fr1", animate=0)
    pymol.cmd.mview("store", 1, 7)
    for f in range(2, num_frames+1):
        pymol.cmd.madd(str(f) + ' x7')
        pymol.cmd.scene("fr" + str(f), animate=0)
        pymol.cmd.mview("store", (f-1)*7+1, f*7)

    pymol.cmd.mview("interpolate")
    pymol.cmd.mview("smooth")
    pymol.cmd.set("movie_loop", 1)
    return
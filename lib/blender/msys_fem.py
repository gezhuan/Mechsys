# Modules

import os
import math
import pickle
import Blender
import bpy
import mechsys as ms
import msys_dict as di

def fill_mesh(obj):
    if obj!=None and obj.type=='Mesh':

        # Mesh
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)

        # MechSys::Mesh::Generic
        mg = ms.mesh_generic()

        # Transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # Vertices
        res = di.get_verts_on_edges_with_tags (obj, msh)
        mg.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            if i in res: onbry = True
            else:        onbry = False
            mg.set_vert (i, onbry, v.co[0], v.co[1])

        # Elements
        nelems = obj.properties['nelems']
        mg.set_nelems (nelems)
        for i in range(nelems):
            mg.set_elem (i, obj.properties['elems']['tags'][i],
                            obj.properties['elems']['onbs'][i],
                            obj.properties['elems']['vtks'][i])
            for j in range(len(obj.properties['elems']['cons'] [str(i)])): mg.set_elem_con  (i, j, obj.properties['elems']['cons'] [str(i)][j])
            for j in range(len(obj.properties['elems']['etags'][str(i)])): mg.set_elem_etag (i, j, obj.properties['elems']['etags'][str(i)][j])
            try:
                for j in range(len(obj.properties['elems']['ftags'][str(i)])): mg.set_elem_ftag (i, j, obj.properties['elems']['ftags'][str(i)][j])
            except:
                pass

        # Restore local coordinates
        msh.verts = ori

        # End
        if edm: Blender.Window.EditMode(1)
        return mg


def run_fea(obj):
    # set cursor
    Blender.Window.WaitCursor(1)

    # mesh
    m = fill_mesh (obj)

    # boundary conditions
    nbrys = di.get_nbrys_numeric (obj)
    ebrys = di.get_ebrys_numeric (obj)
    fbrys = di.get_fbrys_numeric (obj)
    eatts = di.get_eatts_numeric (obj)

    # set geometry
    g = ms.geom(2)
    ms.set_geom (m, nbrys, ebrys, fbrys, eatts, g)

    # solve
    sol = ms.solver('ForwardEuler')
    sol.set_geom(g).set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
    sol.solve()

    # output
    fn = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    ms.out_vtu (g, fn)
    print '[1;34mMechSys[0m: file <'+fn+'> generated'

    # call ParaView
    os.popen('paraview --data='+fn)

    # restore cursor
    Blender.Window.WaitCursor(0)


def gen_script(obj):
    # boundary conditions
    nbrys = di.get_nbrys_numeric (obj)
    ebrys = di.get_ebrys_numeric (obj)
    fbrys = di.get_fbrys_numeric (obj)
    eatts = di.get_eatts_numeric (obj)

    # text
    fn  = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    txt = Blender.Text.New(obj.name+'_script')
    txt.write ('import bpy\n')
    txt.write ('import mechsys\n')
    txt.write ('import msys_fem\n')
    txt.write ('obj   = bpy.data.objects["'+obj.name+'"]\n')
    txt.write ('mesh  = msys_fem.fill_mesh (obj)\n')
    txt.write ('nbrys = '+nbrys.__str__()+'\n')
    txt.write ('ebrys = '+ebrys.__str__()+'\n')
    txt.write ('fbrys = '+fbrys.__str__()+'\n')
    txt.write ('eatts = '+eatts.__str__()+'\n')
    txt.write ('geom  = mechsys.geom(2)\n')
    txt.write ('mechsys.set_geom (mesh, nbrys, ebrys, fbrys, eatts, geom)\n')
    txt.write ('sol   = mechsys.solver("ForwardEuler")\n')
    txt.write ('sol.set_geom(geom)\n')
    txt.write ('sol.set_lin_sol("LA").set_num_div(1).set_delta_time(0.0)\n')
    txt.write ('sol.solve()\n')
    txt.write ('mechsys.out_vtu(geom, "'+fn+'")\n')

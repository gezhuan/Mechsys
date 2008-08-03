# Modules
import math
import time
import pickle
import Blender
import bpy
import mechsys as ms
import msys_dict as di


def print_timing(func):
    def wrapper(*arg):
        t1  = time.time ()
        res = func      (*arg)
        t2  = time.time ()
        print '[1;34mMechSys[0m: %s took [1;31m%f[0m [1;32mseconds[0m' % (func.func_name, (t2-t1))
        return res
    return wrapper


# ======================================================================================= CAD

def add_point(x, y, z):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj==None:
        msh = bpy.data.meshes.new('points')
        obj = scn.objects.new(msh,'points')
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        msh.verts.extend(x,y,z)
        if edm: Blender.Window.EditMode(1)
        obj.select(1)
        Blender.Window.RedrawAll()
    else:
        Blender.Draw.PupMenu('ERROR|Select a Mesh object before calling this function')


def add_points_from_file(filename):
    Blender.Window.WaitCursor(1)
    edm  = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key  = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    file = open(filename, 'r')
    msh  = bpy.data.meshes.new('points')
    scn  = bpy.data.scenes.active
    obj  = scn.objects.new(msh,key)
    for line in file.readlines():
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif words[0]=='x' or words[0]=='X': pass
        else:
            if len(words)>2: x, y, z = float(words[0]), float(words[1]), float(words[2])
            else:            x, y, z = float(words[0]), float(words[1]), 0
            msh.verts.extend(x,y,z)
    if edm: Blender.Window.EditMode(1) # enter edm
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)


def add_spline_from_file(filename):
    Blender.Window.WaitCursor(1)
    edm  = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key  = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    file = open(filename, 'r')
    scn  = bpy.data.scenes.active
    spl  = bpy.data.curves.new(key,'Curve')
    obj  = scn.objects.new(spl,key)
    first = 1
    for line in file.readlines():
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif words[0]=='x' or words[0]=='X': pass
        else:
            if len(words)>2: x, y, z = float(words[0]), float(words[1]), float(words[2])
            else:            x, y, z = float(words[0]), float(words[1]), 0
            point = Blender.BezTriple.New((x, y, z))
            if first==1:
                spl.appendNurb(point)
                path  = spl[0]
                first = 0
            else:
                path.append(point)
    for point in path:
        point.handleTypes = [Blender.BezTriple.HandleTypes.AUTO, Blender.BezTriple.HandleTypes.AUTO]
    spl.update()
    if edm: Blender.Window.EditMode(1) # enter edm
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)


def distance(p1,p2):
    d = p1 - p2
    return math.sqrt(math.pow(d[0],2)+math.pow(d[1],2)+math.pow(d[2],2))


def closest(pt,pointlist):
    idx = 0
    min = distance(pt,pointlist[idx])
    for i in range(1,len(pointlist)):
        d = distance(pt,pointlist[i])
        if d<min:
            min = d
            idx = i
    return min, idx


def arc_point(msh,cen,sp,ep,steps):
    dr1 = sp-cen
    dr2 = ep-cen
    axi = Blender.Mathutils.CrossVecs(dr2,dr1)
    ang = Blender.Mathutils.AngleBetweenVecs(dr2,dr1)
    rot = Blender.Mathutils.RotationMatrix(ang/steps, 3, 'r', axi)
    for i in range(steps-1):
        pt  = cen+rot*dr1
        msh.verts.extend(pt)
        msh.edges.extend(msh.verts[len(msh.verts)-2],msh.verts[len(msh.verts)-1])
        dr1 = pt-cen


def fillet(radius,steps):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        print '[1;34mMechSys:[0m', msh.edges.selected(), 'edges selected' # TODO: check why sometimes 3 edges are selected
        if len(msh.edges.selected())==2:
            # points (vertices)
            v1 = msh.edges[msh.edges.selected()[0]].v1
            v2 = msh.edges[msh.edges.selected()[0]].v2
            v3 = msh.edges[msh.edges.selected()[1]].v1
            v4 = msh.edges[msh.edges.selected()[1]].v2
            # intersection
            i1, i2 = Blender.Mathutils.LineIntersect(v1.co,v2.co,v3.co,v4.co)
            # fillet
            if i1!=i2: Blender.Draw.PupMenu('ERROR|These two edges do not intersect (obj=%s)' % obj.name)
            else:
                #
                #                   ep  _,--* p2
                #                  _,-*'
                #     a=|ep-p1|,--' /
                #        _,--'     /radius  cen
                #    p1 *_--------|----------*-  d=|c-p1|
                #         `--,_    \
                #              `--,_\
                #                   `-*,_
                #                   sp   `--* p3
                di, i = closest(i1,[v1.co,v2.co])
                dj, j = closest(i1,[v3.co,v4.co])
                if i==0:
                    v1.co = i1
                    p1,p2 = v1,v2
                else:
                    v2.co = i1
                    p1,p2 = v2,v1
                if j==0:
                    msh.edges[msh.edges.selected()[1]].v1 = p1
                    if v3!=p1: msh.verts.delete(v3)
                    p3 = msh.edges[msh.edges.selected()[1]].v2
                else:
                    msh.edges[msh.edges.selected()[1]].v2 = p1
                    if v4!=p1: msh.verts.delete(v4)
                    p3 = msh.edges[msh.edges.selected()[1]].v1
                if radius>0.0:
                    # vectors along the edges
                    e1  = p2.co-p1.co; e1.normalize()
                    e2  = p3.co-p1.co; e2.normalize()
                    alp = Blender.Mathutils.AngleBetweenVecs(e1,e2)/2.0
                    a   = radius/math.tan(math.radians(alp))
                    d   = math.sqrt(math.pow(a,2)+math.pow(radius,2))
                    ep  = p1.co+a*e1
                    sp  = p1.co+a*e2
                    mid = Blender.Mathutils.MidpointVecs(sp,ep); mid-=p1.co; mid.normalize()
                    cen = p1.co+d*mid
                    # add points to the arch
                    msh.verts.extend(sp)
                    msh.edges.extend(p3,msh.verts[len(msh.verts)-1])
                    arc_point(msh,cen,sp,ep,steps)
                    msh.verts.extend(ep)
                    msh.edges.extend(msh.verts[len(msh.verts)-2],msh.verts[len(msh.verts)-1])
                    msh.edges.extend(msh.verts[len(msh.verts)-1],p2)
                    msh.verts.delete(p1)
                Blender.Window.RedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select (only) two edges (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1)
    else: Blender.Draw.PupMenu('ERROR|Please, select a Mesh object before calling this function')


def break_edge(at_mid=False):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        if len(msh.edges.selected())==1 or at_mid:
            for e in msh.edges.selected():
                # points (vertices)
                v1 = msh.edges[e].v1
                v2 = msh.edges[e].v2
                # break-point
                if at_mid: bp = 0.5*(v1.co+v2.co)
                else:      bp = Blender.Window.GetCursorPos()
                if bp!=v1.co and bp!=v2.co:
                    msh.verts.extend (bp[0],bp[1],bp[2])
                    new = msh.verts[-1]
                    v2  = msh.edges[e].v2
                    msh.edges[e].v2 = new
                    msh.edges.extend (new,v2)
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1)
    else: Blender.Draw.PupMenu('ERROR|Please, select a Mesh object before calling this function')


# =========================================================================== Structured mesh

#  2+       r_idx    = 0        # start right vertex index
#   |\      edge_ids = [0,1,2]  # indexes of all edges to search for r_idx
#  0| \1    v1_ids   = [1,0,0]  # v1 vertex indexes
#   |  \    v2_ids   = [2,2,1]  # v2 vertex indexes
#  1+---+0  eds      = [1,0,2]  # result: edge indexes
#     2     ids      = [2,1,0]  # result: vertex indexes
def sort_edges_and_verts(msh, edge_ids, r_idx):
    # loop over all connected edges
    eds = []
    ids = []
    while len(edge_ids)>0:
        v1_ids = [msh.edges[ie].v1.index for ie in edge_ids] # ie => index of an edge
        v2_ids = [msh.edges[ie].v2.index for ie in edge_ids]
        if r_idx in v1_ids:
            idx   = v1_ids.index(r_idx)
            r_idx = v2_ids[idx]
            eds.append   (edge_ids[idx])
            ids.append   (r_idx)
            edge_ids.pop (idx)
        elif r_idx in v2_ids:
            idx   = v2_ids.index(r_idx)
            r_idx = v1_ids[idx]
            eds.append   (edge_ids[idx])
            ids.append   (r_idx)
            edge_ids.pop (idx)
        else: break
    return eds, ids


def gen_blk_2d(obj, msh):

    # Local verts indexes:      Global verts idxs:         
    #                                                       
    #     3      6      2          ids[4]   ids[3]   ids[2] 
    #      @-----@-----@                @-----@-----@       
    #      |           |                |           |       
    #      |           |                |           |       
    #    7 @           @ 5       ids[5] @           @ ids[1]
    #      ^y          |                ^y          |       
    #      |  x        |                |  x        |       
    #      @-->--@-----@                @-->--@-----@       
    #    0      4      1          origin   x_plus   ids[0] 
    #                         
    #                         
    # Local edges indexes:      Global edges idxs:
    #                            
    #            y+                    eds[4]   eds[3]
    #      +----(3)----+                @-----@-----@   
    #      |           |                |           |eds[2]
    #      |           |          eds[5]|           |   
    #  x- (0)         (1) x+            @           @
    #      |           |                |           |eds[1]
    #      |           |          eds[6]|           |   
    #      +----(2)----+                @-----@-----@   
    #            y-                    x_axis  eds[0]

    # get local system
    origin, x_plus, y_plus, z_plus = di.get_local_system (obj)
    if origin>-1:
        x_axis = di.get_local_axis (obj, 'x')

        # sort edges
        edge_ids = [e.index for e in msh.edges if e.index!=x_axis]
        eds, ids = sort_edges_and_verts (msh, edge_ids, x_plus)

        # transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # generate lists with coordinates
        ord_ids = [origin, ids[0], ids[2], ids[4], x_plus, ids[1], ids[3], ids[5]]
        cox     = [msh.verts[iv].co[0] for iv in ord_ids]
        coy     = [msh.verts[iv].co[1] for iv in ord_ids]

        # restore local coordinates
        msh.verts = ori

        # generate weights
        nx = di.get_ndiv (obj, 'x')
        ny = di.get_ndiv (obj, 'y')
        if nx<1: nx = 1
        if ny<1: ny = 1
        wx = [1 for i in range(nx)]
        wy = [1 for i in range(ny)]

        # fill lists with tags
        eds   = [x_axis, eds[1], eds[3], eds[5]]
        etags = di.get_tags_list (obj, 'edge', eds)

        # return list with block data
        return [[cox, coy], wx, wy, etags]

    else:
        Blender.Draw.PupMenu('ERROR|Please, define local axes first (obj=%s)' % obj.name)
        return []


def gen_blk_3d(obj, msh):

    # get local system
    origin, x_plus, y_plus, z_plus = di.get_local_system (obj)
    if origin>-1:
        # sort vertices
        vs = []
        ms.block3d_sort (origin, x_plus, y_plus, z_plus, [ed.key for ed in msh.edges], vs)

        # transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # generate list with coordinates
        cox = [msh.verts[iv].co[0] for iv in vs]
        coy = [msh.verts[iv].co[1] for iv in vs]
        coz = [msh.verts[iv].co[2] for iv in vs]

        # restore local coordinates
        msh.verts = ori

        # generate weights
        nx = di.get_ndiv (obj, 'x')
        ny = di.get_ndiv (obj, 'y')
        nz = di.get_ndiv (obj, 'z')
        if nx<1: nx = 1
        if ny<1: ny = 1
        if nz<1: nz = 1
        wx = [1 for i in range(nx)]
        wy = [1 for i in range(ny)]
        wz = [1 for i in range(nz)]

        # return list with block data
        return [[cox, coy, coz], wx, wy, wz]
    else:
        Blender.Draw.PupMenu('ERROR|Please, define local axes first (obj=%s)' % obj.name)
        return []


def gen_struct_mesh():
    # get objects
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)

    # generate blocks
    bks = []
    for obj in obs:
        if obj!=None and obj.type=='Mesh':
            # Mesh
            if len(obj.getAllProperties())==0: Blender.Draw.PupMenu('ERROR|Please, assign all mesh properties to this object(%s) first' % obj.name)
            msh = obj.getData(mesh=1)

            # 2D - linear
            if len(msh.edges)==4:
                print '2D - linear => not yet'

            # 2D - o2
            elif len(msh.edges)==8:
                res = gen_blk_2d (obj, msh)
                if len(res)>0:
                    bks.append    (ms.mesh_block())
                    bks[-1].set2d (di.get_btag(obj), res[0], res[1], res[2])
                    if len(res[3]): bks[-1].set_etags (res[3])
                    obj.select (0)
                else: return

            # 3D - linear
            elif len(msh.edges)==12:
                print '3D - linear => not yet'

            # 3D - o2
            elif len(msh.edges)==24:
                res = gen_blk_3d (obj, msh)
                if len(res)>0:
                    bks.append    (ms.mesh_block())
                    bks[-1].set3d (di.get_btag(obj), res[0], res[1], res[2], res[3])
                    obj.select (0)
                else: return

            else: Blender.Draw.PupMenu('ERROR|Each block must have: 2D=8 edges, 3D=24 edges')

    # generate mesh and draw results
    if len(bks)>0:
        # set cursor
        Blender.Window.WaitCursor(1)

        # generate mesh and write VTU file for ParaView
        mms = ms.mesh_structured()
        ne  = mms.generate (bks)
        fn  = Blender.sys.makename(ext='_MESH_'+obj.name+'.vtu')
        mms.write_vtu (fn)
        print '[1;34mMechSys[0m: [1;33m%d[0m elements generated' % ne
        print '[1;34mMechSys[0m: File <[1;33m%s[0m> created' % fn

        # draw generated mesh
        draw_mesh (mms)

        # restore cursor
        Blender.Window.WaitCursor(0)


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh():
    # get objects
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)

    # set input polygon
    if obj!=None and obj.type=='Mesh':
        Blender.Window.WaitCursor(1)
        msh = obj.getData(mesh=1)
        mu  = ms.mesh_unstructured()
        ets = di.get_tags_ (obj, 'edge')
        rgs = di.get_regs (obj)
        hls = di.get_hols (obj)
        #mu.set_3d(0)
        mu.set_poly_size (len(msh.verts), len(msh.edges), len(rgs), len(hls))
        for v in msh.verts:
            mu.set_poly_point (v.index, v.co[0], v.co[1], v.co[2])
        for e in msh.edges:
            if e.index in ets: mu.set_poly_segment (e.index, e.v1.index, e.v2.index, ets[e.index])
            else:              mu.set_poly_segment (e.index, e.v1.index, e.v2.index)
        print rgs
        print hls
        for i, r in enumerate(rgs):
            mu.set_poly_region (i, int(r[0]), float(r[1]), float(r[2]), float(r[3]))
        for i, h in enumerate(hls):
            mu.set_poly_hole (i, float(h[0]), float(h[1]), float(h[2]))
        maxarea  = di.get_maxarea  (obj)
        minangle = di.get_minangle (obj)
        mu.generate      (float(maxarea), float(minangle))
        draw_mesh (mu)
        Blender.Window.WaitCursor(0)


# ====================================================================================== Draw

@print_timing
def set_elems(obj, nelems, elems):
    obj.properties['nelems'] = nelems
    obj.properties['elems']  = elems


@print_timing
def set_etags(obj, msh, etags):
    for et in etags:
        edge_id = msh.findEdges (et[0], et[1])
        di.set_tag (obj, 'edge', edge_id, etags[et])


@print_timing
def draw_mesh(mms):
    # get vertices and edges
    verts = []
    edges = []
    mms.get_verts (verts)
    mms.get_edges (edges)

    # add new mesh to Blender
    key     = di.get_key()
    scn     = bpy.data.scenes.active
    new_msh = bpy.data.meshes.new      (key+'_structured')
    new_obj = scn.objects.new (new_msh, key+'_structured')
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (1)
    print '[1;34mMechSys[0m: Mesh extended'

    # set elements
    elems  = {}
    nelems = mms.get_elems (elems)
    set_elems (new_obj, nelems, elems)

    # set etags
    etags = {}
    mms.get_etags (etags)
    set_etags (new_obj, new_msh, etags)

    # redraw
    Blender.Window.QRedrawAll()

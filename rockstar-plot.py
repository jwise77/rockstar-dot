import os, sys
import numpy as na
from yt.mods import *

if len(sys.argv) < 3:
    print "usage: %s tree_file halo_id" % sys.argv[0]
    sys.exit()

halo_id = int(sys.argv[-1])
fname = sys.argv[-2]

def load_pfs(fn):
    pf_store = os.path.dirname(fn) + "/../pfs.txt"
    lines = open(pf_store, "r").readlines()
    pfs = {}
    for l in lines:
        if not l.startswith("#"):
            parts = l.split()
            num = int(parts[1])
            f = os.path.dirname(fn) + "/../../" + parts[0]
            pfs[num] = load(f)
    return pfs

def load_redshifts(pfs):
    ka = na.empty(len(pfs), dtype='int')
    za = na.empty(len(pfs), dtype='float64')
    i = 0
    for k,pf in pfs.items():
        ka[i] = k
        za[i] = pf.current_redshift
        i += 1
    isort = ka.argsort()
    result = {'index': ka[isort], 'redshift': za[isort], 'scale': 1.0/(1+za[isort])}
    return result

def load_tree(hid, fn):
    if not os.path.exists(fn):
        print "Cannot find %s" % fn
        sys.exit()
    lines = open(fn, "r").readlines()

    # Search for tree
    start_line = None
    rank = -hid
    num = 0
    if hid <= 0:
        for l in lines:
            if l.startswith("#tree"):
                if num == -hid:
                    hid = int(l.split()[1])
                    break
                num += 1
        print "Halo %d (by mass) has ID %d" % (rank, hid)

    header = "#tree %d\n" % hid
    for i,l in enumerate(lines):
        if l.startswith("#Omega_M"):
            h = float(l.split("=")[3])
        if l == header:
            start_line = i+1
        elif l.startswith("#tree") and start_line != None:
            end_line = i
            break

    # Load data from strings
    nleaves = end_line - start_line
    nfields = len(lines[start_line].split())
    data = na.empty((nleaves, nfields))
    for i in range(start_line, end_line):
        data[i-start_line,:] = map(float, lines[i].split())
    del lines
    return hid, data

def get_halo_prop(data, pf, hid, i):
    return {'hid': hid,
            'pos': data[i,17:20] / pf['mpchcm'],
            'vel': data[i,20:23],
            'vrms': data[i,13],
            'spin': data[i,26],
            'mvir': data[i,9] * pf.hubble_constant,
            'rvir': data[i,11] / pf['kpchcm']}

def find_mmp(data, zz, pfs, hid):
    result = {}
    nleaves = data.shape[0]
    global_mmp = hid
    eps = 1e-3
    for i in range(nleaves):
        mmp = int(data[i,14]) == 1
        hid = int(data[i,1])
        if i > 0: 
            desc = int(data[i,3])
        else:
            desc = hid
        a = data[i,0]
        if desc == global_mmp and mmp:
            global_mmp = hid
            # Find matching scale factor and output
            pos = na.where(na.abs(a - zz['scale'])/a < eps)[0][0]
            pf = pfs[zz['index'][pos]]
            result[zz['index'][pos]] = get_halo_prop(data, pf, hid, i)
    return result

halo_id, data = load_tree(halo_id, fname)
pfs = load_pfs(fname)
zz = load_redshifts(pfs)
mmp = find_mmp(data, zz, pfs, halo_id)

# plot histories
import matplotlib.pyplot as plt
n = len(mmp.keys())
a = na.empty(n)
m = na.empty(n)
for i,k in enumerate(mmp.keys()):
    pos = na.where(zz['index'] == k)[0][0]
    a[i] = zz['scale'][pos]
    m[i] = mmp[pos]['mvir']
plt.plot(a,m)
plt.savefig('HaloMassHistory-%d.png' % halo_id)

# Find max FOV (5*rvir)
fov = 0.0
for k,v in mmp.items():
    fov = max(fov, 5*v['rvir'])

# plot projections
for i,k in enumerate(mmp.keys()):
    pos = na.where(zz['index'] == k)[0][0]
    z = zz['redshift'][pos]
    p = ProjectionPlot(pfs[pos], 'x', 'Density', center = mmp[pos]['pos'],
                       width = (fov, '1'))
    p.annotate_sphere(mmp[pos]['pos'], mmp[pos]['rvir'])
    p.annotate_point(mmp[pos]['pos']+0.51*fov,"z = %.2f" % (1.0/a[i]-1), 
                     text_args = {'size':20})
    p.save('rs-%d-%4.4d' % (halo_id, pos))


import yt
import os, sys, glob
import numpy as np

class consistent_trees():
    def __init__(self, tree_file="./rockstar_halos/trees/tree_0_0_0.dat"):
        if not os.path.exists(tree_file):
            raise RuntimeError("Tree file not found: %s" % (tree_file))
        self.tree_file = tree_file
        self.tree_dir = os.path.dirname(self.tree_file)
        self.rs_dir = os.path.dirname(self.tree_dir)
        self.sim_dir = os.path.dirname(self.rs_dir)
        self.halo_id = None
        lfn = os.path.join(self.tree_dir, "locations.dat")
        # Load tree locations in the file
        if os.path.exists(lfn):
            self.locations = []
            with open(lfn, "r") as fp:
                for l in fp.readlines():
                    if l.startswith("#"): continue
                    values = l.split()
                    self.locations.append([int(values[0]), int(values[2])])
            self.locations = np.array(self.locations)
        else:
            self.locations = None
        # Read header
        header = open(self.tree_file, "r").readline()[1:]
        self.names = header.split()
        # Remove counters in header titles
        for i,n in enumerate(self.names):
            i0 = n.find("(")
            i1 = n.find(")")
            if i0 >= 0 and i1 >= 0:
                instr = n[i0+1:i1]
                if instr.isdigit():
                    self.names[i] = n[:i0]
        # Load the Hubble parameter
        with open(self.tree_file, "r") as fp:
            for l in fp.readlines():
                if l.startswith("#Omega_M"):
                    self.h = float(l.split("=")[3])
                    break
        # Load scale factors of the outputs
        fn = os.path.join(os.path.join(self.rs_dir, "outputs"), "scales.txt")
        self.scale_factors = np.loadtxt(fn)
        self.find_enzo_outputs()
        return

    def find_enzo_outputs(self, bases=None, eps=1e-4):
        # Setup a dict with the rockstar output numbers as the keys
        # and the enzo parameter filenames as the values.
        self.enzo_fn = {}
        _bases = [["DD", "output_"],
                  ["DD", "data"],
                  ["DD", "DD"],
                  ["RD", "RedshiftOutput"],
                  ["RS", "restart"]]
        if bases != None:
            bases = bases + _bases
        else:
            bases = _bases
        all_files = []
        for b in bases:
            all_files += glob.glob("%s/%s????/%s????" % (self.sim_dir, b[0], b[1]))
        for f in all_files:
            ds = yt.load(f)
            scale = 1.0 / (ds.current_redshift+1)
            dela = abs(scale - self.scale_factors[:,1])
            mina = dela.argmin()
            num = int(self.scale_factors[mina,0])
            if dela[mina] < eps:
                self.enzo_fn[num] = f
        return
    
    def line_to_halo(self, line):
        values = map(float, line.split())
        result = {}
        for i in range(len(values)):
            result[self.names[i]] = values[i]
        return result

    def set_halo(self, halo_id):
        """
        Inputs
        ======
        num: If positive number, number of the tree ID.  If negative
        number, select the N-th most massive halo in the latest
        redshift, where N = -num.
        """
        # Search for tree if the locations.dat file doesn't exist.
        if self.locations == None:
            start_line = None
            rank = -halo_id
            num = 0
            if halo_id <= 0:
                with open(self.tree_file, "r") as fp:
                    for l in fp.readlines():
                        if l.startswith("#tree"):
                            if num == -halo_id:
                                halo_id = int(l.split()[1])
                                break
                            num += 1
                print "Halo %d (by mass) has ID %d" % (rank, halo_id)
            self.halo_id = halo_id
            header = "#tree %d\n" % (self.halo_id)
            pos = 0
            with open(self.tree_file, "r") as fp:
                for l in fp.readlines():
                    pos += len(l)
                    if l == header:
                        self.byte_location = pos
                        break
        else:
            if halo_id > 0 and halo_id not in self.locations[:,0]:
                raise RuntimeWarning("halo_id = %d found not in locations.dat" % (halo_id))
                return
            if halo_id <= 0:
                i = -halo_id
                self.halo_id = self.locations[i,0]
            else:
                self.halo_id = halo_id
                i = np.where(self.locations[:,0] == self.halo_id)[0][0]
            self.byte_location = self.locations[i,1]
        return

    def check_halo_set(self):
        if self.halo_id == None:
            raise RuntimeError("Must set halo id with set_halo() first.")
        return

    def get_latest_halo(self):
        self.check_halo_set()
        with open(self.tree_file) as fp:
            fp.seek(self.byte_location)
            line = fp.readline()
        return self.line_to_halo(line)

    def get_all_progenitors(self):
        self.check_halo_set()
        halos = []
        with open(self.tree_file) as fp:
            fp.seek(self.byte_location)
            for l in fp:
                if l.startswith("#tree"): break
                halos.append(self.line_to_halo(l))
        return halos

    def get_enzo_fn(self, scale, eps=1e-5):
        dela = np.abs(scale - self.scale_factors[:,1])
        mina = dela.argmin()
        num = int(self.scale_factors[mina,0])
        return self.enzo_fn[num]
    
    def get_progenitors(self, num=None, z=None, a=None, eps=1e-5):
        self.check_halo_set()
        if num == None and z == None and a == None:
            raise RuntimeError("Must provide at least one of the following: num=, z=, a=\n"
                               "\t num: Number of the output (see scales.txt or self.scale_factor)\n"
                               "\t z:   redshift\n"
                               "\t a:   scale factor")
        if num != None:
            ii = np.where(self.scale_factors[:,0] == num)[0][0]
            if ii == -1:
                raise RuntimeError("Output number %d not found.  See outputs/scales.txt for valid choices")
            a = self.scale_factors[ii,1]
        elif z != None:
            a = 1.0/z - 1
        result = []
        all_progenitors = self.get_all_progenitors()
        for prog in all_progenitors:
            dela = prog['scale'] - a
            if abs(dela) < eps:
                result.append(prog)
        return result

    def get_mm_progenitors(self):
        self.check_halo_set()
        # Go through the progenitors and select the branches that
        # contain the halo's most-massive progenitor in the complete
        # tree, starting with the latest output.
        mm_chain = []
        for num in range(int(self.scale_factors[-1,0]),
                         int(self.scale_factors[0,0])-1, -1):
            prog = self.get_progenitors(num=num)
            for p in prog:
                if bool(p['mmp?']):
                    mm_chain.append(p)
                    break
        return mm_chain

    def list_scale_factors(self):
        for a in self.scale_factors:
            print "Output %4d: a = %.6f" % (a[0], a[1])
        return

if __name__ == "__main__":
    fn = sys.argv[-1]
    if not os.path.exists(fn):
        raise RuntimeError("Tree file not found: %s" % (fn))
    ct = consistent_trees(fn)
    ct.set_halo(0)
    mm = ct.get_mm_progenitors()
    for p in mm:
        print "(z = %.2f) Mvir = %.4g Msun" % (1/p['scale']-1, p['mvir']/ct.h)

from __future__ import print_function
from pylab import *
from fsl.caches import *
import nibabel as nb


class MelodicIOError(IOError):
    pass


def trace_call(fn):
    from itertools import chain
    import traceback

    def wrapped(*v, **k):
        name = fn.__name__
        print("{0}({1})".format(name, ", ".join(map(repr, chain(v, k.values())))))
        traceback.print_stack(limit=2)
        return fn(*v, **k)
    return wrapped


class melodic:
    def __init__(self, path="", ic=1):
        self.set_directory(path)
        self.select_component(ic)

    def set_directory(self, path=""):
        self.dir = path
        self.ic = 1

        mfname = "{0}/melodic_mix".format(self.dir)
        ffname = "{0}/melodic_FTmix".format(self.dir)
        bfname = "{0}/mean.nii.gz".format(self.dir)
        vfname = "{0}/melodic_ICstats".format(self.dir)
        zfname = "{0}/melodic_IC.nii.gz".format(self.dir)

        self.mix = genfromtxt(mfname)
        self.FTmix = genfromtxt(ffname)
        self.exp_var_stats = genfromtxt(vfname)
        self.background_image = bfname
        self.zstat_data = nb.load(zfname).get_data()
        self.select_component(self.ic)

    def select_component(self, ic=1):
        self.ic = ic

        self.stat = self.zstat_data[:, :, :, self.ic-1]

    def get_tr(self):
        tr = 1
        for line in open("{0}/log.txt".format(self.dir)):
            if "--tr" in line:
                print(line)
                for arg in line.split():
                    if "--tr" in arg:
                        tr = float(arg.split("=")[1])
                break
        return tr

    def get_bg_image(self):
        return self.background_image

    def get_mix_shape(self):
        return self.mix.shape

    def get_stat_data(self):
        return self.stat

    def get_mix_data(self):
        return self.mix[:, self.ic-1]

    def get_FTmix_data(self):
        return self.FTmix[:, self.ic-1]

    def get_variance_stats(self):
        return (self.exp_var_stats[self.ic-1, 0], self.exp_var_stats[self.ic-1, 1])

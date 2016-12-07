from __future__ import print_function
import nibabel as nb
import numpy as np
from os.path import join


class Melodic:
    def __init__(self, path="", ic=1):
        self._set_directory(path)
        self.select_component(ic)

    def _set_directory(self, path):
        self.dir = path
        self.ic = 1

        mfname = join(self.dir, "melodic_mix")
        ffname = join(self.dir, "melodic_FTmix")
        bfname = join(self.dir, "mean.nii.gz")
        #vfname = join(self.dir, "melodic_ICstats")
        zfname = join(self.dir, "melodic_IC.nii.gz")

        self.mix = np.genfromtxt(mfname)
        self.FTmix = np.genfromtxt(ffname)
        #self.exp_var_stats = np.genfromtxt(vfname)
        self.background_image = bfname
        self.zstat_data = nb.load(zfname).get_data()
        self.select_component(self.ic)

    def select_component(self, ic):
        self.ic = ic
        self.stat = self.zstat_data[:, :, :, self.ic - 1]

    @property
    def tr(self):
        trval = 1
        for line in open(join(self.dir, "log.txt")):
            if "--tr" in line:
                print(line)
                for arg in line.split():
                    if "--tr" in arg:
                        trval = float(arg.split("=")[1])
                break
        return trval

    @property
    def bg_image(self):
        return self.background_image

    @property
    def mix_shape(self):
        return self.mix.shape

    @property
    def stat_data(self):
        return self.stat

    @property
    def mix_data(self):
        return self.mix[:, self.ic - 1]

    @property
    def FTmix_data(self):
        return self.FTmix[:, self.ic - 1]

    #def get_variance_stats(self):
    #    return (self.exp_var_stats[self.ic-1, 0], self.exp_var_stats[self.ic-1, 1])

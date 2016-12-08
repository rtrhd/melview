#!/usr/bin/env python

"""
    Melview

    Copyright(c) 2012, University of Oxford (David Flitney)

    This file is part of Melview

    Melview is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division, print_function, absolute_import
from os.path import join, split, expanduser, exists, dirname, isabs, isdir
from os import makedirs, environ
import re
import csv

environ['ETS_TOOLKIT'] = 'qt4'
environ['QT_API'] = 'pyside'

import logging
logging.basicConfig(level=logging.INFO)

import numpy as np
import nibabel as nb
from nibabel.py3k import FileNotFoundError

from . melodic import Melodic
try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt

from pyface.api import FileDialog, OK, confirm, error, YES

from matplotlib.figure import Figure
from matplotlib import cm

from traits.api import (
    HasTraits, File, Directory, List,
    Button, Str, Int, Float, Bool, Enum, Instance, Range,
    on_trait_change
)
from traitsui.api import (
    Handler, View, Item,
    HGroup, HSplit, VSplit, VGroup, Group,
    RangeEditor
)


# etsconfig moved to traits since ETS 4.0.
# If enthought is not installed try to import etsconfig from traits.
try:
    from enthought.etsconfig.api import ETSConfig
except ImportError:
    from traits.etsconfig.api import ETSConfig

from traitsui.menu import Action
from traitsui.api import TableEditor
from traitsui.table_column import ObjectColumn


if ETSConfig.toolkit == "qt4":
    from traitsui.qt4.editor import Editor
    from traitsui.qt4.basic_editor_factory import BasicEditorFactory

from mpl_toolkits.axes_grid1 import make_axes_locatable

from nibabel.spatialimages import ImageFileError

from configobj import ConfigObj
import validate
import platform
from errno import EEXIST


if platform.system() == 'Linux':
    settings_file = join(expanduser('~'), '.config', 'melview', 'melviewrc')
else:
    settings_file = join(expanduser('~'), '.melviewrc')
if not exists(dirname(settings_file)):
    try:
        makedirs(dirname(settings_file))
    except OSError as e:  # Guard against race condition
        if e.errno != EEXIST:
            raise

opts = [
    "plot_height = float(0, 2400, default=400)",
    "plot_width = float(0, 1, default=0.8)",
    "height = float(0, 1, default=0.8)",
    "width = float(0, 1, default=0.8)",
    "ncols = integer(0, 1000, default=0)"
]
settings = ConfigObj(settings_file, configspec=opts)

validator = validate.Validator()
settings.validate(validator, copy=True)


class _MPLFigureEditor(Editor):

    scrollable = True

    def init(self, parent):
        self.control = self._create_canvas(parent)

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        if ETSConfig.toolkit == "qt4":
            from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
            from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as _NavigationToolbar
            from matplotlib.backends.backend_qt4 import cursord
            from pyface.qt import QtGui

            class NavigationToolbar(_NavigationToolbar):
                def set_cursor(self, cursor):
                    QtGui.QApplication.restoreOverrideCursor()
                    qcursor = QtGui.QCursor()
                    qcursor.setShape(cursord[cursor])
                    QtGui.QApplication.setOverrideCursor(qcursor)

            # with additional widgets such as the Matplotlib navigation toolbar
            widget = QtGui.QWidget()
            widget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)

            canvas = FigureCanvas(self.value)
            canvas.setParent(widget)
            toolbar = NavigationToolbar(canvas, widget)
            vbox = QtGui.QVBoxLayout(widget)
            vbox.addWidget(canvas)
            vbox.addWidget(toolbar)

            return widget


class MPLFigureEditor(BasicEditorFactory):
    klass = _MPLFigureEditor


class Classification(HasTraits):
    class_names = [
        'Signal', 'Unknown', 'Unclassified Noise', 'Movement',
        'Cardiac', 'White Matter', 'Non-brain', 'MRI', 'Susceptability-motion',
        'Sagittal sinus', 'Respiratory'
    ]
    filter      = Bool(False)
    dirty       = Bool(False)
    class_name  = Enum(values='class_names')
    ic_number   = Int
    display_min = Float(0.0)
    display_max = Float(0.0)

    def __init__(self, **kwargs):
        super(Classification, self).__init__(**kwargs)
        self.dirty = False

    def __str__(self):
        return "{0}, {1}, {2}".format(self.ic_number, self.class_name, self.filter)

    @on_trait_change('filter, class_name')
    def changed(self, b):
        self.dirty = True
        self.filter = self.class_name not in ('Signal', 'Unknown')

    view = View(Item('class_name', show_label=False, style='custom'))


class ClassificationColumn(ObjectColumn):

    def get_cell_color(self, object):
        return 'light blue' if object.class_name == 'Signal' else 'light green'


class MelodicWindowHandler(Handler):

    def do_quit(self, info):
        info.ui.dispose()

    def close(self, info, isok):
        if info.object.dirty():
            response = confirm(info.ui.control, """
Results have not been saved yet.

Are you sure you still want to exit?
""")
            return response == YES
        else:
            return True

    def do_load(self, info):
        if info.object.dirty():
            response = confirm(info.ui.control, """
Loading will overwrite unsaved results.

Are you sure you still want to continue?
""")

        if not info.object.dirty() or (response == YES):
            dialog = FileDialog(action="open", wildcard=info.object.file_wildcard)
            dialog.open()
            if dialog.return_code == OK:
                info.object.filedir = dialog.directory
                info.object.filename = dialog.filename
                info.object.file_load()

            del dialog

    def do_save(self, info):
        info.object.file_save()

    def do_save_as(self, info):
        default_directory = split(info.object.dirpath)[0]
        dialog = FileDialog(
            action="save as", wildcard=info.object.file_wildcard,
            default_directory=default_directory
        )
        dialog.open()
        if dialog.return_code == OK:
            info.object.filedir = dialog.directory
            info.object.filename = dialog.filename
            info.object.file_save()

        del dialog


class MelodicWindow(HasTraits):

    dirpath = Directory
    bgimage = File
    initialised = False
    threshold = Float(2.3)

    filepath = ""
    filedir = ""

    class_list = List(Classification)

    lbox = Instance(Figure, ())
    plots = Instance(Figure, ())

    nics = Int(1)
    ic = Int(0)

    tr = Float(2.0)

    lut_min = Float(0)
    lut_max = Float(10)
    background_min = Float(0)
    background_max = Float(10)

    ncols = Int(settings['ncols'])

    needs_saving = Bool(False)

    show_stats = Bool(True)
    ignore_blank_slices = Bool(True)

    filedir = Str
    filename = Str
    file_wildcard = Str("Text file (*.txt)|*.txt|Data file (*.dat)|*.dat|All files|*")

#    quit    = Action(name="Quit", action='do_quit')

    load = Action(name="Load classifications...", action="do_load")
    save = Action(
        name="Save classifications", action="do_save",
        enabled_when='needs_saving and filename is not ""'
    )
    save_as = Action(name="Save as...", action="do_save_as", enabled_when='needs_saving')

    reset_button = Button(label="Reset LUT to data max/min")
    record_zoom = Button(label="Record Zoom")
    reset_zoom = Button(label="Reset Zoom")

    from traitsui.key_bindings import KeyBinding, KeyBindings

    key_bindings = KeyBindings(
        KeyBinding(
            binding1='Ctrl-A',
            description='previous IC',
            method_name='_prev_button_fired'
        ),
        KeyBinding(
            binding1='Ctrl-Z',
            description='next IC',
            method_name='_next_button_fired'
        )
    )

    ic_selected = Instance(Classification)

    ic_list_editor = TableEditor(
        columns=[
            ClassificationColumn(name='ic_number', label='IC#', editable=False),
            ClassificationColumn(name='class_name')
        ],
        selected='ic_selected',
        selection_mode='row',
        rows=5,
        row_factory=Classification
    )

    orientation_names = ['XY', 'XZ', 'YZ']
    view_orientation = Enum(values='orientation_names')

    view = View(
        Group(
            VSplit(
                HSplit(
                    Item(
                        'lbox',
                        editor=MPLFigureEditor(),
                        show_label=False, width=settings['plot_width'], height=settings['plot_height']
                    ),
                    Item(
                        'class_list',
                        editor=ic_list_editor,
                        show_label=False, springy=False, resizable=False
                    )
                ),
                HSplit(
                    Item(
                        'plots',
                        editor=MPLFigureEditor(),
                        show_label=False, width=0.7, height=0.2
                    ),
                    VGroup(
                        Item('ic_selected', label='Classification', style='custom'),
                        HGroup(
                            Item('show_stats', label='Show statistics'),
                            Item('lut_min', label="LUT min/max"),
                            Item('lut_max', show_label=False)
                        ),
                        Item('ignore_blank_slices'),
                        HGroup(Item('reset_button', show_label=False)),
                        HGroup(
                            Item('record_zoom', show_label=False),
                            Item('reset_zoom', show_label=False),
                            Item('ncols', editor=RangeEditor(Range(0, 100), mode='spinner'))
                        ),
                        HGroup(Item('view_orientation'))
                    )
                )
            ),
            label='Explorer'
        ),
        Group(
            Item('dirpath', label='ICA directory'),
            Item('bgimage', label='Structural image'),
            Item('tr'),
            Item('threshold'),
            Item('background_min'),
            Item('background_max'),
            label='Data'
        ),
        title='Melodic Results Viewer',
        width=settings['width'],
        height=settings['height'],
        resizable=True,
        key_bindings=key_bindings,
        buttons=[load, save, save_as],
        handler=MelodicWindowHandler()
    )

    def montageXZ(self, X, ncols=0):
        x, y, z = X.shape
        if ncols == 0:
            mm = int(np.ceil(np.sqrt(z) * 1.5) - 1)
        else:
            mm = ncols
        nn = y // mm + 1

        M = np.zeros((mm * x, nn * z))

        image_id = 0
        for j in range(mm):
            for k in range(nn):
                if image_id >= y:
                    break
                sliceM, sliceN = j * x, k * z
                M[sliceM:sliceM + x, sliceN:sliceN + z] = X[:, image_id, :]
                image_id += 1

        return M

    def montageYZ(self, X, ncols=0):
        x, y, z = X.shape
        if ncols == 0:
            mm = int(np.ceil(np.sqrt(y) * 1.5) - 1)
        else:
            mm = ncols
        nn = x // mm + 1

        M = np.zeros((mm * y, nn * z))

        image_id = 0
        for j in range(mm):
            for k in range(nn):
                if image_id >= x:
                    break
                sliceM, sliceN = j * y, k * z
                M[sliceM:sliceM + y, sliceN:sliceN + z] = X[image_id, :, :]
                image_id += 1

        return M

    def montageXY(self, X, ncols=0):
        x, y, z = X.shape
        if ncols == 0:
            mm = int(np.ceil(np.sqrt(x) * 1.5) - 1)
        else:
            mm = ncols
        nn = z // mm + 1
        M = np.zeros((mm * x, nn * y))

        image_id = 0
        for j in range(mm):
            for k in range(nn):
                if image_id >= z:
                    break
                sliceM, sliceN = j * x, k * y
                M[sliceM:sliceM + x, sliceN:sliceN + y] = X[:, :, image_id]
                image_id += 1

        return M

    def __init__(self, path="", bgpath="", ic=1, *args, **kwargs):
        logging.debug("MelodicWindow.__init__()")
        super(MelodicWindow, self).__init__(*args, **kwargs)

        self.montage = self.montageXY
        try:
            self.image_axes = self.lbox.add_subplot(1, 1, 1)
            self.image_axes.hold(True)
            self.lbox.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)

            self.mix_axes = self.plots.add_subplot(1, 2, 1)
            self.mix_axes.hold(False)

            self.ftmix_axes = self.plots.add_subplot(1, 2, 2)
            self.ftmix_axes.hold(False)
            self.plots.subplots_adjust(right=0.95, bottom=0.2)

            self.bg_data = nb.load(bgpath).get_data()
            self.dirpath = path
            self.bgimage = bgpath

            self.mel.select_component(1)
            if self.initialised:
                self.display()
        except (FileNotFoundError, ImageFileError, AttributeError):
            pass
        finally:
            self.initialised = True

        logging.debug("MelodicWindow.__init__()...done")

    def dirty(self):
        unsaved = False
        for f in self.class_list:
            if f.dirty:
                unsaved = True
                break
        self.needs_saving = unsaved
        return unsaved

    def mark_clean(self):
        for f in self.class_list:
            f.dirty = False
        self.needs_saving = False

    @lru_cache(maxsize=10)
    def get_bgdata(self, path):
        return nb.load(path).get_data()

    def nonzero_slice_list(self, X):
        nx, ny, nz = X.shape
        return [n for n in range(nz) if np.max(X[:, :, n]) > 0]

    @on_trait_change('bgimage')
    def bgimage_changed(self, obj, name, oldpath, path):
        logging.debug("name={0}, oldpath={1}, path={2}".format(name, oldpath, path))
        try:
            self.bg_data = self.get_bgdata(path)
            self.background_min, self.background_max = 0, self.bg_data.max()
        except IOError as e:
            error(
                None,
                """
Couldn't change the background image!

Attempting to revert to previous value.

{0}
""".format(e),
                title="Melodic file error"
            )
            if oldpath:
                self.bgimage = oldpath
        if self.initialised:
            self.display()

    @on_trait_change('dirpath')
    def dirpath_changed(self, obj, name, oldpath, path):
        try:
            self.mel = Melodic(path, ic=1)
            self.bgimage = self.mel.bg_image
            self.update_info()
            self.tr = self.mel.tr
            if self.initialised:
                self.display()

            cl = []
            for ic in range(self.nics):
                cl.append(
                    Classification(
                        ic_number=ic+1,
                        class_name='Unknown',
                        filter=False
                    )
                )

            self.class_list[:] = cl

            self.reset_lut()
        except IOError as e:
            error(
                None, """
    Couldn't find all the required data in {0}.

    Perhaps this isn't a Melodic directory?

    {1}
    """.format(path, e),
                title="Melodic dir error"
            )

    def reset_lut(self):
        ic = self.ic_selected
        if ic is not None:
            if ic.display_max == ic.display_min == 0.0:
                self.lut_min, self.lut_max = self.threshold, self.mel.stat_data.max()
            else:
                self.lut_min, self.lut_max = ic.display_min, ic.display_max

    def update_info(self):
        npts, nics = self.mel.mix_shape
        self.nics = nics

    @on_trait_change('lut_min,lut_max')
    def change_lut(self):
        ic = self.ic_selected
        if ic is not None and self.initialised:
            if self.lut_max > self.lut_min:
                ic.display_min, ic.display_max = self.lut_min, self.lut_max
                self.display()

    @on_trait_change('background_min,background_max')
    def change_bglut(self):
        if self.initialised:
            if self.background_max > self.background_min:
                self.display()

    def display(self):
        logging.debug("MelodicWindow.display()")
        d = self.mel.stat_data
        b = self.bg_data

        if self.ignore_blank_slices:
            slices = self.nonzero_slice_list(b)

            zi = np.rot90(self.montage(d[:, :, slices], ncols=self.ncols))
            bg = np.rot90(self.montage(b[:, :, slices], ncols=self.ncols))
        else:
            zi = np.rot90(self.montage(d, ncols=self.ncols))
            bg = np.rot90(self.montage(b, ncols=self.ncols))

        self.image_axes.clear()
        self.image_axes.axis('off')
        bgim = self.image_axes.imshow(
            bg, cmap='gray',
            vmin=self.background_min,
            vmax=self.background_max
        )
        bgim.set_interpolation("nearest")

        if self.show_stats and self.ic_selected is not None:
            cm1 = cm.hot
            cm1.set_under(alpha=0)
            im1 = self.image_axes.imshow(zi, cmap=cm1, vmin=self.lut_min, vmax=self.lut_max)
            im1.set_interpolation("nearest")
            cm2 = cm.Blues
            cm2.set_over(alpha=0)
            im2 = self.image_axes.imshow(zi, cmap=cm2, vmax=-self.lut_min, vmin=-self.lut_max)
            im2.set_interpolation("nearest")

        if self.ic_selected is not None:
            if not hasattr(self, "colorbar"):
                divider = make_axes_locatable(self.image_axes)
                cax = divider.append_axes("right", size="2.5%", pad=0.05)
                self.colorbar = self.image_axes.get_figure().colorbar(im1, cax=cax)
            self.colorbar.set_clim(self.lut_min, self.lut_max)
            self.colorbar.update_normal(d)

            if self.ic_selected.class_name == 'Signal':
                self.image_axes.set_title(
                    64*' ' + '%d [Signal]' % self.ic_selected.ic_number + 80*' ',
                    bbox={'color': 'lightblue'}, fontsize=20
                )
            else:
                self.image_axes.set_title(
                    '%d [%s]' % (self.ic_selected.ic_number, self.ic_selected.class_name),
                    fontsize=20
                )

        y = self.mel.mix_data
        t = self.tr * np.arange(len(y))
        self.mix_axes.plot(t, y)
        self.mix_axes.set_xlabel('Time (s)')

        y = self.mel.FTmix_data
        n = len(y)
        f = 1.0 / (2 * n * self.tr) * np.arange(n)
        self.ftmix_axes.plot(f, y)
        self.ftmix_axes.set_xlabel('Frequency (Hz)')

        if self.ftmix_axes.figure.canvas is not None:
            self.ftmix_axes.figure.canvas.draw()

        if self.initialised:
            self.lbox.canvas.draw()

        logging.debug("MelodicWindow.display()... done")

    @on_trait_change('ncols')
    def ncolsChange(self, ncols):
        self.display()

    @on_trait_change('view_orientation')
    def orientationChange(self, orient):
        self.montage = getattr(self, 'montage{0}'.format(orient))
        self.display()

    @on_trait_change('threshold')
    def thresholdChanged(self, threshold):
        self.lut_min = threshold
        self.display()

    @on_trait_change('tr')
    def trChanged(self, tr):
        if tr > 0.0:
            self.display()

    @on_trait_change('class_list')
    def classListChanged(self):
        self.display()

    @on_trait_change('show_stats,ignore_blank_slices')
    def optionsChanged(self):
        self.image_axes.clear()
        self.display()

    @on_trait_change('ic_list_editor')
    def selectedElementChanged(self):
        pass
        #self.display()

    @on_trait_change('ic_selected')
    def selectedChanged(self, obj):
        limits = self.image_axes.axis()
        if obj:
            self.mel.select_component(obj.ic_number)

        self.needs_saving = self.dirty()
        self.reset_lut()
        self.reset_axes(limits)

    def reset_axes(self, limits):
        self.image_axes.axis(limits)
        self.image_axes.figure.canvas.draw()

    def _record_zoom_fired(self, info):
        self.alimits = self.image_axes.axis()

    def _reset_zoom_fired(self, info):
        if hasattr(self, 'alimits'):
            self.reset_axes(self.alimits)

    def _prev_button_fired(self, info):
        if self.ic_selected.ic_number > 1:
            ic = self.ic_selected.ic_number - 2
            self.ic_selected = self.class_list[ic]

    def _reset_button_fired(self, info):
        self.lut_min, self.lut_max = self.threshold, self.mel.stat_data.max()
        self.reset_lut()

    def _next_button_fired(self, info):
        if self.ic_selected.ic_number < self.nics:
            ic = self.ic_selected.ic_number
            self.ic_selected = self.class_list[ic]

    def file_load(self):
        logging.debug("MelodicWindow.file_load()")
        path = join(self.filedir, self.filename)
        logging.debug("Loading file {0}".format(path))

        try:
            with open(path) as f:
                lines = [line.strip() for line in f]
            if not len(lines) < 2:
                raise ValueError('Specified file (%s) has less than two lines' % path)

            icadirpath = lines[0] if isabs(lines[0]) else join(self.filedir, lines[0])
            if not isdir(icadirpath):
                raise ValueError('Directory path (%s) specified in %s does not exist' % (icadirpath, path))

            components = [
                Classification(ic_number=ic+1, class_name='Unknown', filter=False)
                for ic in range(self.nics)
            ]

            rows = [(int(i), c, f == 'True') for i, c, f in csv.reader(lines[1:-1])]
            if not all(len(row) == 3 for row in rows):
                raise ValueError('Badly formed component file "%s"' % path)
            for i, c, f in rows:
                components[i-1] = Classification(ic_number=i, class_name=c, filter=f)

            do_fix = len(rows) < 1
            if do_fix:
                last_line = lines[-1]
                if re.match('^\[[0-9, ]*\]$', last_line):
                    for icstr in last_line[1:-1].split(','):
                        components[int(icstr)-1].class_name = 'Unclassified Noise'
                else:
                    raise ValueError("Bad input in last line of classifications file!")

            self.dirpath = icadirpath
            self.class_list[:] = components
            self.ic_selected = self.class_list[0]

            self.mark_clean()
            self.display()

        except (IOError, ValueError) as e:
            error(None, 'Unable to read classifications file {0}\n{1}'.format(path, e), title='File Error')

        logging.debug("MelodicWindow.file_load()...done")

    def file_save(self):
        path = join(self.filedir, self.filename)
        try:
            with open(path, 'w') as f:
                print(self.dirpath, file=f)
                for c in self.class_list:
                    print(c, file=f)
                print('[' + ','.join([int(ic.ic_number) for ic in self.class_list if ic.filter]) + ']', file=f)
        except IOError as e:
            error(None, 'Unable to open file for writing\n{0}\n{1}'.format(path, e), title='File I/O error')

        self.mark_clean()
        self.display()


def main():
    cfg = [
        "path   = string(default='')",
        "bgpath = string(default='')"
    ]

    config = ConfigObj('melview.ini', configspec=cfg)
    config.validate(validator, copy=True)
    m = MelodicWindow(path=config['path'], bgpath=config['bgpath'])
    m.configure_traits()

    config['path'] = m.dirpath
    config['bgpath'] = m.bgimage

    try:
        config.write()
    except (FileNotFoundError, IOError) as e:
        error(None, 'Error while writing config file to working directory (%s).' % e, title='Config error')

    settings['ncols'] = m.ncols
    try:
        settings.write()
    except (FileNotFoundError, IOError) as e:
        error(None, 'Error while writing settings file to home directory (%s).' % e, title='Settings error')


if __name__ == "__main__":
    main()

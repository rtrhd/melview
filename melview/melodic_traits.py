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

import os
import re
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyside'

import logging
logging.basicConfig(level=logging.INFO)

import numpy as np

from fsl.melodic import *

#import Qt4

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt
#from pyface.qt import QtGui, QtCore

from pyface.api import FileDialog, OK, confirm, error, YES

#from PyQt4.QtGui import QWidget, QSizePolicy

import matplotlib
from matplotlib.figure import Figure 

from traits.api import \
    HasTraits, File, Directory, List, \
    Button, Str, Int, Float, Bool, Enum, Any, Instance, Range,\
    on_trait_change
from traitsui.api import  Handler, View, Item, \
    UItem, HGroup, HSplit, VSplit, VGroup, Group, \
    RangeEditor, EnumEditor, InstanceEditor

from traitsui.api import CustomEditor
from enthought.etsconfig.api import ETSConfig

from traitsui.menu import Action
if ETSConfig.toolkit == "wx":
    from traitsui.wx.editor import Editor
    from traitsui.wx.basic_editor_factory import BasicEditorFactory
elif ETSConfig.toolkit == "qt4":
    from traitsui.qt4.editor import Editor
    from traitsui.qt4.basic_editor_factory import BasicEditorFactory

from mpl_toolkits.axes_grid1 import make_axes_locatable

from nibabel.spatialimages import ImageFileError
#import gc

from configobj import ConfigObj
import validate

validator = validate.Validator()
opts= """
plot_height = float(default=400)
plot_width = float(default=0.8)
height = float(default=0.8)
width = float(default=0.8)
ncols = integer(default=0)
"""
settings_file = os.path.join(os.path.expanduser('~'), '.melviewrc')
settings = ConfigObj(settings_file, configspec = opts.split("\n"))
settings.validate(validator, copy=True)

class _MPLFigureEditor(Editor):

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
#        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        if ETSConfig.toolkit == "wx":
            from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
            from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
            import wx
                    
            """ Create the MPL canvas. """
            # The panel lets us add additional controls.
            panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
            sizer = wx.BoxSizer(wx.VERTICAL)
            panel.SetSizer(sizer)
            # matplotlib commands to create a canvas
            mpl_control = FigureCanvas(panel, -1, self.value)
            sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
            toolbar = NavigationToolbar(mpl_control)
            sizer.Add(toolbar, 0, wx.EXPAND)
            self.value.canvas.SetMinSize((10,10))
        
            return panel

        elif ETSConfig.toolkit == "qt4":
            from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
            from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as _NavigationToolbar
            from matplotlib.backends.backend_qt4 import cursord
            from pyface.qt import QtGui
            
            class NavigationToolbar(_NavigationToolbar):
                def set_cursor( self, cursor ):
                    QtGui.QApplication.restoreOverrideCursor()
                    qcursor = QtGui.QCursor()
                    qcursor.setShape(cursord[cursor])
                    QtGui.QApplication.setOverrideCursor( qcursor )
        
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


#def montage(X, colormap=cm.gist_gray):    
#    m, n, count = shape(X)
#    mm = int(ceil(sqrt(count) * 1.5)-1)
#    nn = count / mm + 1
##     print m, n, count
##     print mm, nn
#    M = zeros((mm * m, nn * n))
#
#    image_id = 0
#    for j in range(mm):
#        for k in range(nn):
#            if image_id >= count: 
#                break
#            sliceM, sliceN = j * m, k * n
##  	    print j, k, sliceM, sliceN
#            M[sliceM:sliceM + m, sliceN:sliceN + n] = X[:, :, image_id]
#            image_id += 1
#                    
#    return M

class Classification(HasTraits):
    class_names= ['Signal', 'Unknown', 'Unclassified Noise', 'Movement', \
		  'Cardiac', 'White Matter', 'Non-brain', 'MRI', 'Susceptability-motion', \
                  'Sagittal sinus','Respiratory']
    filter     = Bool(False)
    dirty      = Bool(False)
    class_name = Enum(values='class_names')
    ic_number  = Int
    display_min = Float(0.0)
    display_max = Float(0.0)

    def __init__(self, **kwargs):
        super( Classification, self).__init__( **kwargs )
        self.dirty = False

    def __str__(self):
        return "{0}, {1}, {2}".format(self.ic_number, self.class_name, self.filter)

    @on_trait_change('filter, class_name')
    def changed(self, b):
        self.dirty = True
        self.filter = self.class_name not in ('Signal', 'Unknown')

    view = View( Item('class_name', show_label=False, style='custom') )

from traitsui.api import TableEditor
from traitsui.table_column import ObjectColumn
from enthought.traits.api import Color
from enthought.traits.ui.api import TabularEditor

class ClassificationColumn(ObjectColumn):
    
    def get_cell_color ( self, object ):
        return [ 'light blue', 'light green' ][ object.class_name in ['Signal']]

#class CheckboxColumn(ObjectColumn):
#    
#    def __init__ ( self, **traits ):
#        """ Initializes the object.
#        """
#        super( CheckboxColumn, self ).__init__( **traits )
#    #---------------------------------------------------------------------------
#    #  Returns the cell background color for the column for a specified object: 
#    #---------------------------------------------------------------------------
#   
#    def get_cell_color ( self, object ):
#        """ Returns the cell background color for the column for a specified
#            object.
#        """
#        # we override this from the parent class to ALWAYS provide the
#        # standard color
#        return self.cell_color_
#    #---------------------------------------------------------------------------
#    #  Returns whether the column is editable for a specified object: 
#    #---------------------------------------------------------------------------
#               
#    def is_editable ( self, object ):
#        """ Returns whether the column is editable for a specified object.
#        """
#        # Although a checkbox column is always editable, we return this
#        # to keep a standard editor from appearing. The editing is handled
#        # in the renderer's handlers.
#        return False

class MelodicWindowHandler(Handler):
    
    def do_quit(self,info):
#        self.closesession()
        info.ui.dispose()

    def close(self, info, isok):
        if info.object.dirty():
            response = confirm(info.ui.control,
"""
Results have not been saved.

Are you sure you want to exit?
""")
            return response == YES
        else:
            return True

    def do_load(self, info):
        if info.object.dirty():
            response = confirm(info.ui.control,
"""
Loading will overwrite unsaved results.

Are you sure you want to continue?
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

#    def do_reset_lut(self, info):
#        info.object.reset_lut()
#
#    def do_record_zoom(self, info):
#        info.object.reset_zoom()
#        
#    def do_reset_zoom(self, info):
#        info.object.reset_zoom()
        
    def do_save_as(self, info):
        default_directory = os.path.split(info.object.dirpath)[0]
        dialog = FileDialog(action="save as", wildcard=info.object.file_wildcard,
                            default_directory=default_directory)
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
    filedir  = ""

    class_list = List(Classification)

    lbox   = Instance(Figure, ())
    plots  = Instance(Figure, ())

    nics = Int(1)
    ic   = Int(0)

    tr   = Float(2.0)

    lut_min = Float(0)
    lut_max = Float(10)
    background_min = Float(0)
    background_max = Float(10)

    ncols = Int(settings['ncols'])
       
    needs_saving = Bool(False)

    show_stats = Bool(True)
    ignore_blank_slices = Bool(True)

    filedir  = Str
    filename = Str
    file_wildcard = Str("Text file (*.txt)|*.txt|Data file (*.dat)|*.dat|All files|*")

#    quit    = Action(name="Quit", action='do_quit')

    load    = Action(name = "Load classifications...", action = "do_load")
    save    = Action(name = "Save classifications", action = "do_save",
             enabled_when='needs_saving and filename is not ""')
    save_as = Action(name = "Save as...", action = "do_save_as", enabled_when='needs_saving')

#    reset = Action(name="Reset LUT", action='do_reset_lut')

    reset_button = Button(label="Reset LUT to data max/min")
    record_zoom  = Button(label="Record Zoom")
    reset_zoom   = Button(label="Reset Zoom")
    
    from traitsui.key_bindings import KeyBinding, KeyBindings

    key_bindings = KeyBindings(
    KeyBinding( binding1 = 'Ctrl-A',
            description = 'previous IC',
            method_name = '_prev_button_fired' ),
    KeyBinding( binding1 = 'Ctrl-Z',
            description = 'next IC',
            method_name = '_next_button_fired' )
    )

    ic_selected = Instance(Classification)

    ic_list_editor = TableEditor(columns = [ ClassificationColumn( name='ic_number', label='IC#', editable=False ), 
                                             ClassificationColumn( name='class_name' ) ],
                                selected         = 'ic_selected',
                                selection_mode   = 'row',
                                rows             = 5,
                                row_factory      = Classification)
    
    orientation_names = ['XY', 'XZ', 'YZ']
    view_orientation = Enum(values='orientation_names')

    view = View(
                Group(
                      VSplit(
                             HSplit(
                                    Item('lbox',
                                         editor=MPLFigureEditor(),
                                         show_label=False, width=settings['plot_width'], height=settings['plot_height']),
                                    Item('class_list',
                                         editor=ic_list_editor,
                                         show_label=False, springy=False, resizable=False)
                                    ),
                             HSplit(
                                    Item('plots',
                                         editor=MPLFigureEditor(),
                                         show_label=False, width=0.7, height=0.2),
                                    VGroup(Item('ic_selected', label='Classification', style='custom'),
                                           HGroup(Item('show_stats', label='Show statistics'),
                                                  Item('lut_min', label="LUT min/max"),
                                                  Item('lut_max', show_label=False)),
                                           Item('ignore_blank_slices'),
                                           HGroup(Item('reset_button', show_label=False)),
                                           HGroup(Item('record_zoom', show_label=False),
                                                  Item('reset_zoom', show_label=False),
                                                  Item('ncols', editor=RangeEditor(Range(0, 100), mode='spinner'))),
                                           HGroup(Item('view_orientation'))
                                           )
                                    ) 
                             ),
                      label='Explorer'),
                Group(Item('dirpath', label='ICA directory'), 
                      Item('bgimage', label='Structural image'),
                      Item('tr'),
                      Item('threshold'),
                      Item('background_min'),
                      Item('background_max'),
                      label='Data'),
                title   = 'Melodic Results Viewer',
                width   = settings['width'],
                height  = settings['height'],
#                x = 100, y = 100,
                resizable    = True,
                key_bindings = key_bindings,
#        buttons = [ load, save, save_as, quit ],
                buttons = [ load, save, save_as ],
                handler = MelodicWindowHandler()
                )

    def montageXZ(self, X, ncols=0):
        x, y, z = shape(X)
        if ncols == 0:
            mm = int(ceil(sqrt(z) * 1.5) - 1)
        else:
            mm = ncols
        nn = y / mm + 1
        
        M = zeros((mm * x, nn * z))
    
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
        x, y, z = shape(X)
        if ncols == 0:
            mm = int(ceil(sqrt(y) * 1.5) - 1)
        else:
            mm = ncols
        nn = x / mm + 1
        
        M = zeros((mm * y, nn * z))
    
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
        x, y, z = shape(X)
        if ncols == 0:
            mm = int(ceil(sqrt(x) * 1.5) - 1)
        else:
            mm = ncols
        nn = z / mm + 1
        M = zeros((mm * x, nn * y))

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
#        print "MelodicWindow.__init__() called super"
                
        self.montage = self.montageXY
        try:
            self.image_axes = self.lbox.add_subplot(1,1,1)
            self.image_axes.hold(True)
            self.lbox.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)   
#            print "MelodicWindow.__init__() created image_axes"

            self.mix_axes = self.plots.add_subplot(1,2,1)
            self.mix_axes.hold(False)
#            print "MelodicWindow.__init__() created mix_axes"

            self.ftmix_axes = self.plots.add_subplot(1,2,2)
            self.ftmix_axes.hold(False)
            self.plots.subplots_adjust(right=0.95, bottom=0.2)   
#            print "MelodicWindow.__init__() created ftmix_axes"

#            print "MelodicWindow.__init__() reading data"
            self.bg_data = nb.load(bgpath).get_data()
            self.dirpath = path
            self.bgimage = bgpath
#            print "MelodicWindow.__init__() read data"

            self.mel.select_component(1)
            if self.initialised:
                self.display()
        except ImageFileError, e:
            pass
        except Exception, e:
            print e
            raise e
        finally:
            self.initialised = True

        logging.debug("MelodicWindow.__init__()...done")

    def dirty(self):
        unsaved = False
        for f in self.class_list :
            if f.dirty:
                unsaved = True
                break
        self.needs_saving = unsaved
        return unsaved

    def mark_clean(self):
        for f in self.class_list :
            f.dirty = False
        self.needs_saving = False

    @lru_cache(maxsize=10)
    def get_bgdata(self, path):
        return nb.load(path).get_data()
    
    def nonzero_slice_list(self, X):
        (nx, ny, nz) = X.shape
        return [n for n in range(0, nz) if np.max(X[:,:,n]) > 0]

    @on_trait_change('bgimage')
    def bgimage_changed(self, obj, name, oldpath, path):
        print "name={0}, oldpath={1}, path={2}".format(name, oldpath, path)
        try:
            self.bg_data = self.get_bgdata(path) 
            self.background_min, self.background_max = 0, self.bg_data.max()
        except IOError, e:
            error(None, 
    """
    Couldn't change the background image!
    
    Attempting to revert to previous value.
    
    {0}
    """.format(e), 
    		  title="Melodic file error")
            if oldpath != "":
                self.bgimage = oldpath
        if self.initialised:
            self.display()

    @on_trait_change('dirpath')
    def dirpath_changed(self, obj, name, oldpath, path):
        try:
            self.mel = melodic(path, ic=1)
            self.bgimage = self.mel.get_bg_image()
            self.update_info()
            self.tr = self.mel.get_tr()
            if self.initialised:
                self.display()
    
            cl = []
            for ic in arange(0, self.nics):
                cl.append( Classification(ic_number=ic+1, 
                                          class_name='Unknown', 
                                          filter=False) )
    
            self.class_list[:] = cl
    
            self.reset_lut()
        except IOError, e:
            error(None, 
    """
    Couldn't find all the required data in {0}. 
    
    Perhaps this isn't a Melodic directory?
    
    {1}
    """.format(path, e),
    		  title="Melodic dir error")
#             self.dirpath=oldpath

    def reset_lut(self):
        ic = self.ic_selected
        if ic:
            if ic.display_max == 0.0 and ic.display_min == 0.0:
                self.lut_min, self.lut_max = self.threshold, self.mel.get_stat_data().max()
            else:
                self.lut_min, self.lut_max = ic.display_min, ic.display_max
    
    def update_info(self):
        npts, nics = self.mel.get_mix_shape()
    #	print "update_info nics = {0}".format(nics)
        self.nics = nics

    @on_trait_change('lut_min,lut_max')
    def change_lut(self):
    #	self.colorbar.set_clim(self.lut_min, self.lut_max)
    #	print "{0} {1}".format(self.lut_min, self.lut_max)
        ic = self.ic_selected
        if self.initialised:
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
        d=self.mel.get_stat_data()
        b=self.bg_data

        if self.ignore_blank_slices:
            slices = self.nonzero_slice_list(b)
        
            zi = rot90(self.montage(d[:,:,slices],ncols=self.ncols))
            bg = rot90(self.montage(b[:,:,slices],ncols=self.ncols))
        else:
            zi = rot90(self.montage(d,ncols=self.ncols))
            bg = rot90(self.montage(b,ncols=self.ncols))
    
        self.image_axes.clear()
        self.image_axes.axis('off')
        bgim = self.image_axes.imshow(bg, cmap=cm.gray, 
    				      vmin=self.background_min, 
    				      vmax=self.background_max)
        bgim.set_interpolation("nearest")
        
        if self.show_stats and not (self.ic_selected == None):
            cm1 = cm.hot
            cm1.set_under(alpha=0)
            im1 = self.image_axes.imshow(zi, cmap=cm1, vmin=self.lut_min, vmax=self.lut_max)
            im1.set_interpolation("nearest")
            cm2 = cm.Blues
            cm2.set_over(alpha=0)
            im2 = self.image_axes.imshow(zi, cmap=cm2, vmax=-self.lut_min, vmin=-self.lut_max)
            im2.set_interpolation("nearest")
    
        if not self.ic_selected == None:
            if not hasattr(self, "colorbar"):
                divider = make_axes_locatable(self.image_axes)
                cax = divider.append_axes("right", size="2.5%", pad=0.05)
                self.colorbar = self.image_axes.get_figure().colorbar(im1, cax=cax)
            self.colorbar.set_clim(self.lut_min, self.lut_max)
            self.colorbar.update_normal(d)

            if self.ic_selected.class_name in 'Signal':
                self.image_axes.set_title('                                                                  Signal                                                                               ', bbox={'color':'lightgreen'})
            else:
                self.image_axes.set_title(self.ic_selected.class_name)
    
        y = self.mel.get_mix_data()
        n = y.shape[0]
        t = arange(0.0, n * self.tr, self.tr)
        self.mix_axes.plot(t[:n], y)
    
        y = self.mel.get_FTmix_data()
        n = y.shape[0]
        f = 1.0 / (2 * n * self.tr)
        t = arange(0.0, n * f, f)
        self.ftmix_axes.plot(t[:n], y)
        if self.ftmix_axes.figure.canvas != None:
            self.ftmix_axes.figure.canvas.draw()
        
        if self.initialised:
            self.lbox.canvas.draw()
            
        logging.debug("MelodicWindow.display()... done")

    @on_trait_change('ncols')
    def ncolsChange(self, ncols):
        self.display()
        
    @on_trait_change('view_orientation')
    def orientationChange(self, orient):
        self.montage=getattr(self, 'montage{0}'.format(orient))
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
    #	print "Ding!"
        self.display()
        pass

    @on_trait_change('show_stats,ignore_blank_slices')
    def optionsChanged(self):
        self.image_axes.clear()
        self.display()
        
#    @on_trait_change('ignore_blank_slices')
#    def ignoreBlanksChanged(self):
#    #    print "Ding!"
#        self.image_axes.clear()
#        self.display()
#    
#    @on_trait_change('show_stats')
#    def showStatsChanged(self):
#    #	print "Ding!"
#        self.image_axes.clear()
#        self.display()

    @on_trait_change('ic_list_editor')
    def selectedElementChanged(self):
    #	print "Ding!"
        pass

    @on_trait_change('ic_selected')
    def selectedChanged(self, obj):
    # 	print "{0}".format(obj)
        limits = self.image_axes.axis()
        
        if obj:
            self.mel.select_component(obj.ic_number)
#            if self.initialised:
#                self.display()
                
        self.needs_saving = self.dirty()
        self.reset_lut()
        self.reset_axes(limits)

    def reset_axes(self, limits):
#        print "reset_axes"
        self.image_axes.axis(limits)
        self.image_axes.figure.canvas.draw()
        
    def _record_zoom_fired(self, info):
        self.alimits = self.image_axes.axis()

    def _reset_zoom_fired(self, info):
        if hasattr(self,'alimits'):
            self.reset_axes(self.alimits)

    def _prev_button_fired(self, info):
    #	print "Prev"
        if self.ic_selected.ic_number > 1:
            ic = self.ic_selected.ic_number - 2
            self.ic_selected = self.class_list[ic]
    
    def _reset_button_fired(self, info):
        self.lut_min, self.lut_max = self.threshold, self.mel.get_stat_data().max()
        self.reset_lut()

    def _next_button_fired(self, info):
    #	print "Next"
        if self.ic_selected.ic_number < self.nics:
            ic = self.ic_selected.ic_number
            self.ic_selected = self.class_list[ic]

    def file_load(self):
        print "MelodicWindow.file_load()"
        path = os.path.join(self.filedir, self.filename)
        print "Loading file {0}".format(path)
    
        try:
            import re
            f = open(path, 'r')
    
            cl = []
    #	    print len(cl)
            lines = f.readlines()
    #	    print lines[0]
    
            icadirpath=lines[0].strip()
            if re.match('^/', icadirpath):
                self.dirpath = icadirpath
            else:
                self.dirpath = os.path.join(self.filedir, icadirpath)

            for ic in arange(0, self.nics):
                cl.append( Classification(ic_number=ic+1, class_name='Unknown', filter=False) )

            doFix = True
            for i,fv,c in [ (int(i),(fv.strip() == 'True'),c.strip()) 
                           for (i,c,fv) in [l.split(',') for l in lines[1:-1]] ] :
    #		print i, fv, c
                cl[i-1] = Classification(ic_number=i, class_name=c, filter=fv)
                doFix = False

            last_line=lines[-1]
            
            if doFix:
                if re.match('\[[0-9, ]*\]\n$', last_line):
                    
                    for ic in eval(last_line):
                        cl[ic-1].class_name='Unclassified Noise'
    
                else:
                    raise Exception("Bad input in last line of classifications file!")
                        
            self.class_list[:] = cl

            self.ic_selected = self.class_list[0]
    
            self.mark_clean()
            self.display()
    
        except e:
            print e
            error(None, 'Unable to open file for reading\n{0}'.format(path), title='File I/O error')

        print "MelodicWindow.file_load()...done"
    
    def file_save(self):
        path = os.path.join(self.filedir, self.filename)
    
        try:
            f = open(path, 'w')
    
            f.write(self.dirpath + '\n')
#            for c in [ c for c in self.class_list if c.class_name != 'Unknown' ]:
            for c in self.class_list:
                f.write( '{0}\n'.format(c) )	
            f.write( str([ int(ic.ic_number) for ic in self.class_list if ic.filter == True]) + '\n')
    
            f.close()

            self.mark_clean()
            self.display()
    
        except e:
            print e
            error(None, 'Unable to open file for writing\n{0}'.format(path), title='File I/O error')
    
def main():
 
#    print "__main__"
    cfg = """
    path   = string(default='')
    bgpath = string(default='')
    """

#    print "__main__ reading config"
    config = ConfigObj('melview.ini', configspec=cfg.split("\n"))
    config.validate(validator, copy=True)
    
#    print "__main__ read config"

    m = MelodicWindow(path=config['path'], bgpath=config['bgpath'])
    m.configure_traits()

    config['path']   = m.dirpath
    config['bgpath'] = m.bgimage

    try:
        config.write()
    except:
        error(None, 'Error while writing config file to working directory.', title='Config error')

    settings['ncols'] = m.ncols
    
    try:
        settings.write()
    except:
        error(None, 'Error while writing settings file to home directory.', title='Settings error')

        
if __name__ == "__main__":
    main()

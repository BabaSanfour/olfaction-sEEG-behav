# -*- coding: utf-8 -*-


from OpenElectrophy.gui.contextmenu import  MenuItem, Delete
from OpenElectrophy.gui.editdb import EditFieldsDialog

from summarywidget import SummaryWidget
from mainviewer import MainViewer


from PyQt4.QtCore import *
from PyQt4.QtGui import *

class Delete(MenuItem):
    name = 'Delete'
    table = None
    mode = 'homogeneous'
    icon = ':/user-trash.png'
    def execute(self, session, treeview, ids, tablename,  treedescription, **kargs):
        for warn in  [  'Do you want to delete this and all of its descendants?',
                                'Are you sure?',
                                'Are you really sure?',
                                ]:
            mb = QMessageBox.warning(treeview,treeview.tr('delete'),treeview.tr(warn), 
                    QMessageBox.Ok ,QMessageBox.Cancel  | QMessageBox.Default  | QMessageBox.Escape,
                    QMessageBox.NoButton)
            if mb == QMessageBox.Cancel : return
        
        for id in ids:
            class_ = treedescription.tablename_to_class[tablename]
            instance = session.query(class_).get(id)
            session.delete(instance)
        session.commit()
        treeview.refresh()

class Edit(MenuItem):
    name = 'Edit'
    table = None
    mode = 'unique'
    icon = ':/view-form.png'
    def execute(self, session,treeview, id, tablename,treedescription, explorer,  **kargs):
        class_ = treedescription.tablename_to_class[tablename]
        instance = session.query(class_).get(id)
        w= EditFieldsDialog(parent = treeview, session = session, instance = instance)
        if w.exec_():
            explorer.refresh()

class ViewRunSummary(MenuItem):
    name = 'View Run Summary'
    table = 'Run'
    mode = 'unique'
    icon = ':/plot.png'
    def execute(self, session, id, tablename, treedescription,treeview,  **kargs):
        class_ = treedescription.tablename_to_class[tablename]
        instance = session.query(class_).get(id)
        w = SummaryWidget(parent = treeview, run = instance)
        w.setWindowFlags(Qt.Window)
        w.show()


from OpenElectrophy.gui.viewers import SegmentViewer, SignalViewer, TimeFreqViewer

class DrawSegment(MenuItem):
    name = 'View Segment'
    table = 'Segment'
    mode = 'unique'
    icon = ':/draw-segment.png'
    def execute(self, session,treeview, id, tablename,treedescription, explorer,  **kargs):
        class_ = treedescription.tablename_to_class[tablename]
        neo_seg = session.query(class_).get(id).to_neo(cascade = True)
        w= SegmentViewer(segment = neo_seg, parent = treeview)
        w.setWindowFlags(Qt.Window)
        w.show()


class SignalViewer(MenuItem):
    name = 'Un viewer contre du vin'
    table = 'Run'
    mode = 'unique'
    icon = ':/draw-segment.png'
    def execute(self, session,treeview, id, tablename,treedescription, explorer,  **kargs):
        class_ = treedescription.tablename_to_class[tablename]
        run = session.query(class_).get(id)
        w= MainViewer(run = run, parent = treeview)
        w.setWindowFlags(Qt.Window)
        w.show()

context_menu = [ Delete , ViewRunSummary, Edit, DrawSegment, SignalViewer]

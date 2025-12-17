
# -*- coding: utf-8 -*-


from  connection import *


from PyQt4.QtCore import *
from PyQt4.QtGui import *


from OpenElectrophy import *
from OpenElectrophy.gui import *
from OpenElectrophy.gui.qtsqltreeview import *
from OpenElectrophy.gui.guiutil import *


from contextmenu import context_menu
from insertdata import insert_one_run

import sqlalchemy as sa

treedescription = TreeDescription(dbinfo = dbinfo,
                                    table_children = { 
                                                            'Subject' : ['Run', 'Block',   ],
                                                            'Block' : [ 'RecordingChannelGroup',  ],
                                                            'RecordingChannelGroup' : ['RecordingChannel', ],
                                                            #~ 'RecordingChannel' : ['AnalogSignal'],
                                                            #~ 'Run' : [ 'Trial', ],
                                                            'Run' : [ 'Trial',  'Segment',  ], #'epocharray', 'eventarray', ],
                                                             'Segment' : [ 'AnalogSignal', 'EpochArray', 'EventArray', 'RespirationSignal'],
                                                            },
                                    columns_to_show = { },
                                    table_on_top = 'Subject',
                                    table_order = { 'Run' : asc(sa.text('`Run`.`date`')),
                                                            } ,
                                    )



class MainWindow(QMainWindow) :
    def __init__(self, parent = None,):
        QMainWindow.__init__(self, parent)

        self.resize(900, 600)
        
        self.mainWidget = QWidget()
        self.setCentralWidget(self.mainWidget)
        self.mainLayout = QHBoxLayout()
        self.mainWidget.setLayout(self.mainLayout)

        self.setWindowTitle('Epidosic analyses')
        self.setWindowIcon(QIcon(':/openelectrophy.png'))
        
        
        self.treeview = QtSqlTreeView(session = session, treedescription = treedescription, context_menu = context_menu)
        self.mainLayout.addWidget(self.treeview)
        
        self.createActions()
        self.createMenus()
        
        
        
    def createActions(self):
        self.actionAddFile = QAction('Add file (.res)', self)
        self.actionAddFile.setIcon(QIcon(':/document-open-folder.png'))
        self.actionAddFile.triggered.connect(self.insertData)
        

    def createMenus(self):
        self.fileMenu = self.menuBar().addMenu("File")
        self.fileMenu.addAction(self.actionAddFile)
        #~ self.fileMenu.addSeparator()
        
        #~ self.fileMenu.addAction(self.quitAct)
    
    def insertData(self):
        class Parameters(DataSet):
            filenames = FilesOpenItem("Open files", "res", '')
            channel_trig_odor = StringItem("channel_trig_odor", '')
        p = Parameters()
        if p.edit():
            for filename_res in p.filenames:
                print filename_res
                try:
                    #~ insert_one_run(filename_res)
                    insert_one_run(filename_res, channel_trig_odor = p.channel_trig_odor)
                except:
                    mb = QMessageBox.warning(self,self.tr("Oups"), self.tr("y'a un truc qui cloche"), 
                            QMessageBox.Ok ,QMessageBox.Cancel  | QMessageBox.Default  | QMessageBox.Escape,
                            QMessageBox.NoButton)
                    
            self.treeview.refresh()



def startMain():
    app = QApplication([ ])
    w = MainWindow()
    w.show()
    app.exec_()



if __name__ =='__main__':
    startMain()




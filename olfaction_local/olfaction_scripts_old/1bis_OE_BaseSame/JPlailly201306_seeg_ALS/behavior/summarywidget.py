# -*- coding: utf-8 -*-

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from connection import *


from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar2
from matplotlib.figure import Figure

#~ from respiration_features import compute_cycle_for_run
from respiration_tools import compute_cycle_for_run


class SimpleCanvas(FigureCanvas):
    def __init__(self, parent=None, ):
        self.fig = Figure()
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                    QSizePolicy.Expanding,
                    QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        color = self.palette().color(QPalette.Background).getRgb()
        color = [ c/255. for c in color[:3] ]
        self.fig.set_facecolor(color)




class SimpleCanvasAndTool(QWidget):
    def __init__(self  , parent = None , ):
        QWidget.__init__(self, parent)
        self.mainLayout = QVBoxLayout()
        self.setLayout(self.mainLayout)

        self.canvas = SimpleCanvas()
        self.fig = self.canvas.fig
        self.toolbar = NavigationToolbar2(self.canvas , parent = self ,  )
        
        self.mainLayout.addWidget(self.toolbar)
        self.mainLayout.addWidget(self.canvas)
        


class SummaryWidget(QWidget):
    def __init__(self  , parent = None ,
                            run = None):
        QWidget.__init__(self, parent)
        self.run = run
        
        self.mainLayout = QVBoxLayout()
        self.setLayout(self.mainLayout)
        
        t = 'Summary for subject={} exp={} run={} group={}'.format(run.subject.name, run.exp, run.index, run.group) 
        self.setWindowTitle(t)
        self.setWindowIcon(QIcon(':/openelectrophy.png'))
        
        self.mainLayout.addWidget(QLabel(t))

        self.canvas = SimpleCanvasAndTool()
        self.mainLayout.addWidget(self.canvas)
        
        self.table = QTableWidget()
        self.mainLayout.addWidget(self.table)
        
        self.table.itemSelectionChanged.connect(self.zoom_on_trial)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        
        
        
        self.refresh()

    def refresh(self):
        resp = self.run.respirationsignals[0]
        #odors = self.run.epocharrays[0] #FIXME
        odors = self.run.epocharrays.filter_by(name = 'odors')[0]

        
        times = resp.t_start + np.arange(resp.signal.shape[0])/resp.sampling_rate
        self.canvas.fig.clear()
        ax = self.canvas.fig.add_subplot(111)
        self.ax = ax
        ax.plot(times, resp.signal, color = 'b')
        
        self.resp_f = medfilt(resp.signal, kernel_size = medfilt_kernel_size)
        #~ print self.resp_f
        ax.plot(times, self.resp_f,  color = 'c')
        
        ax.axhline(mean(self.resp_f), color = 'm', linewidth = 2)
        
        
        # compute and plot cycles
        
        if resp.cycle_times is None:
            'compute resp'
            resp.cycle_times = compute_cycle_for_run(resp)
            resp.update(session)
        cycle_times = resp.cycle_times.magnitude

        
        cycles_ind = zeros(cycle_times.shape, dtype = 'i')
        for i in range(2):
            cycles_ind[:, i] = digitize(cycle_times[:,i], times)
        
        ax.plot( cycle_times[:,0], resp.signal[cycles_ind[:,0] ],ls = 'None', marker = 'o', color = 'r')
        ax.plot( cycle_times[:-1,1],resp.signal[cycles_ind[:-1,1] ], ls = 'None', marker = 'o', color = 'g')
        
        
        
        
        
        ax.axhline(resp.odor_threshold, color ='g')
        
        trials = self.run.trials.order_by(asc(text('`Trial`.`index`'))).all()
        
        
        for t1, t2 in zip(odors.times, odors.durations):
            ax.axvspan(t1,t1+t2, color ='r', alpha = 0.3)
        
        for trial in trials:
            ax.axvline(trial.time, color ='k')
            ax.axvline(trial.wanted_odor_time, color = 'r')
            
            od = trial.odor_name#.decode('ascii', 'ignore')
            #~ print od
            ax.text(trial.triggered_odor_time, 5, u'{}'.format(od))
            
            if trial.recognition_time is not None:
                ax.axvline(trial.time + trial.recognition_time, color = 'c')
                ax.text(trial.time + trial.recognition_time, 5, trial.recognition)
            
            if trial.recognition == 1:
                ax.axvline(trial.time + trial.context_time, color = 'g')
                
                ax.axvline(trial.time + trial.click_image_time, color = 'g')
                ax.text(trial.time + trial.click_image_time, 5, str(trial.selected_image))
                ax.axvline(trial.time + trial.click_posistion_time, color = 'g')
                ax.text(trial.time + trial.click_posistion_time, 5, str(trial.click_posistion))

        attrs = Trial.usable_attributes
        
        self.table.setColumnCount(len(attrs))
        self.table.setRowCount(len(trials))
        self.table.setHorizontalHeaderLabels(attrs.keys())
        for row, trial in enumerate(trials):
            for col, attr in enumerate(attrs):
                #~ print attr, (u'{}'.format(getattr(trial, attr))).decode('ascii', 'ignore')
                self.table.setItem(row, col,  QTableWidgetItem(u'{}'.format(getattr(trial, attr))))
            
        
    def zoom_on_trial(self):
        i= self.table.selectedIndexes()[0].row()
        
        
        trials = self.run.trials.order_by(asc(text('`Trial`.`index`'))).all()
        t1 = trials[i].time-5
        if i <len(trials)-1:
            t2 = trials[i+1].time+2
        else:
            t2 = trials[i].time+40
        self.ax.set_xlim(t1,t2)
        
        self.canvas.canvas.draw()
    


def test_SummaryWidget():
    run = session.query(Run).filter_by(exp='E').first()
    #~ run = session.query(Run).get(2)
    app = QApplication([ ])
    w = SummaryWidget(run = run)

    w.show()
    app.exec_()



if __name__ =='__main__':
    test_SummaryWidget()


    

        
    

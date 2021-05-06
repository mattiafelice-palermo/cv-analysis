#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 09:42:10 2021

@author: mattia
"""

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        self.setStyleSheet("background-color:transparent;")
        self.figure.set_facecolor("none")
        self.figure.tight_layout(pad = 0.1)
        
        self.animation = False
       
        self.mouse_pos = [[0, 0], [0, 0]]
        
    def static_plot(self, x, y, labels = ['x', 'y']):
        if self.animation is True:
            dir(self._ani.event_source)
            self._ani.event_source.stop()
            self.animation = False
        
        self.axes.cla()
        self.axes.set_xlabel(labels[0])
        self.axes.set_ylabel(labels[1])
        self.axes.plot(x,y)
        self.figure.tight_layout(pad = 0.1)
        self.draw()        
        
    def animated_plot(self, x, y):
        self.plots = {"base": self.axes.plot([], [])[0],
                      "mouse": self.axes.scatter([], [])
                      }
        self._ani = animation.FuncAnimation(
            self.figure, self._animate, init_func=self._init_frame, interval=40, blit=True
        )
        self.figure.canvas.mpl_connect("motion_notify_event", self._on_mouse_move)
        
        
    def _init_frame(self):
        self.plots['mouse'].set_offsets([0, 0])
        return (self.plots['mouse'],)
    
    def _animate(self, i):
        self._ani.event_source.stop()
        self.axes.relim()
        self.axes.autoscale_view()
        self._ani.event_source.start()
        self._mouse_pointer()
        return (self.plots['mouse'], self.plots['base'])
    
    
    def _on_mouse_move(self, event):
        # element 0 is the most recent mouse pos, element 1 is the previous
        if event.xdata is None:
            self.mouse_pos[0] = self.mouse_pos[1]
        else:
            self.mouse_pos[1] = self.mouse_pos[0]
            self.mouse_pos[0] = [event.xdata, event.ydata]
            
    def _mouse_pointer(self):
        mouse = self.mouse_pos[0][0], self.mouse_pos[0][1]
        self.plots["mouse"].set_offsets([mouse[0], mouse[1]])
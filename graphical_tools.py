#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 15:45:19 2021

@author: mattia
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from sklearn.metrics import r2_score

from scipy.optimize import curve_fit

#%% RangeSelectBase
class RangeSelectBase:
    """
    Base class for creating interactive animations with matplotlib on 2D plots.

    RangeSelectBase provides some basic tools such as selected range highlight,
    mouse pointer, line follower and mouseclick catching.

    Attributes
    ----------
    data : numpy.ndarray
        Stacked numpy array of input x,y data
    fig : matplotlib.figure.Figure
        Figure object from matplotlib
    axes : dict
        Contains minimum and maximum values of the axes (keys: xmin, xmax,ymin, ymax)
    ax : matplotlib.axes._subplots.AxesSubplot
        Matplotlib axes object, must be used to plot data.
    mouse_pos : list
        Coords of current (mouse_pos[0]) and previous (mouse_pos[1]) mouse position
    clicks: list
        list of coordinates of mouse clicks
    plot_open : bool
        True if plot is open, otherwise false
    range_mask : list
        List of bools, True for data within range selected by user
    self.plots : list
        Can contain different artist objects (lines.Line2D, collections.PathCollection, etc)

    Methods
    -------
    run()
        Method description
    custom_animation()
        Method description
    return_function()
        Method description
    run_condition()
        Method description
    """

    def __init__(self, x, y):
        # Data to be plotted and figure
        self.data = np.vstack((x, y)).T
        self.fig = plt.figure()
        self._ani = None

        # define axes limits and create axis
        self.axes = self._set_axes()
        self.ax = plt.axes(
            xlim=(self.axes["xmin"], self.axes["xmax"]),
            ylim=(self.axes["ymin"], self.axes["ymax"]),
        )

        # interactive variables from matplotlib events
        self.mouse_pos = [[0, 0], [0, 0]]
        self.clicks = []
        self.plot_open = True

        # Array containing True/False mask for the selected range
        self.range_mask = None

        # dict of base plots. If they are not fed any data, they just won't display
        self.plots = {
            "base": self.ax.plot([], [])[0],
            "mouse": self.ax.scatter([], []),
            "line follow": self.ax.scatter([], []),
            "range hlight": self.ax.scatter([], []),
        }

    def _set_axes(self, c=0.1):
        xmin, xmax = np.amin(self.data[:, 0]), np.amax(self.data[:, 0])
        ymin, ymax = np.amin(self.data[:, 1]), np.amax(self.data[:, 1])

        xvar = (xmax - xmin) * c
        yvar = (ymax - ymin) * c

        return {
            "xmin": xmin - xvar,
            "xmax": xmax + xvar,
            "ymin": ymin - yvar,
            "ymax": ymax + yvar,
        }

    def _init_frame(self):
        self.plots["base"].set_data(self.data[:, 0], self.data[:, 1])
        return (self.plots["base"],)

    def _animate(self, i):
        self.custom_animation()

        return tuple(self.plots.values())  # tuple of plot objects

    def _mouse_pointer(self):
        mouse = self.mouse_pos[0][0], self.mouse_pos[0][1]
        self.plots["mouse"].set_offsets([mouse[0], mouse[1]])

    def _line_follow(self):
        ref = np.array(self.mouse_pos[0])
        distance = np.linalg.norm(ref - self.data, axis=1)
        index = distance.argmin()
        point = self.data[index]

        self.plots["line follow"].set_offsets([point[0], point[1]])

        return point

    def _range_highlight(self):
        data = self.data
        clicks = self.clicks

        point = self._line_follow()

        if len(clicks) == 1:
            click_dist = data[:, 0] - clicks[0][0]
            point_dist = data[:, 0] - point[0]
            if point[0] < clicks[0][0]:
                mask = np.logical_and(click_dist < 0, point_dist > 0)
            else:
                mask = np.logical_and(click_dist > 0, point_dist < 0)

            self.range_mask = mask
            self.plots["range hlight"].set_offsets(data[mask])

    def _on_mouse_move(self, event):
        # element 0 is the most recent mouse pos, element 1 is the previous
        if event.xdata is None:
            self.mouse_pos[0] = self.mouse_pos[1]
        else:
            self.mouse_pos[1] = self.mouse_pos[0]
            self.mouse_pos[0] = [event.xdata, event.ydata]

    def _on_mouse_lclick(self, event):
        # usability while zooming - to be improved with ESC catch?
        if plt.get_current_fig_manager().toolbar.mode != "":
            return
        self.clicks.append([event.xdata, event.ydata])

    def _on_close(self, event):
        self.plot_open = False

    def run(self):
        # Run graphical analysis process
        self._ani = animation.FuncAnimation(
            self.fig, self._animate, init_func=self._init_frame, interval=40, blit=True
        )

        # connection to events, send signal to functions
        plt.connect("motion_notify_event", self._on_mouse_move)
        plt.connect("button_press_event", self._on_mouse_lclick)
        self.fig.canvas.mpl_connect("close_event", self._on_close)

        # plt.ion() # not really necessary?
        plt.show()

        # run until user configurated condition is False
        while self.run_condition():
            plt.pause(0.1)

        # Return data from the return function!
        return self.return_function()

    def custom_animation(self):
        pass

    def return_function(self):
        pass

    def run_condition(self):
        pass


#%%
class MaxMinFinder(RangeSelectBase):
    def __init__(self, x, y, mode):
        super().__init__(x, y)
        self.plots["show maxima"] = self.ax.scatter([], [], color="red")
        self.maxima = None
        self.max_index = None
        self.mode = mode

    def custom_animation(self):
        self._mouse_pointer()
        self._line_follow()
        self._range_highlight()
        mask = self.range_mask

        data = self.data

        if self.mode == "max":
            function = np.argmax
        if self.mode == "min":
            function = np.argmin

        if len(self.clicks) == 1:
            if len(data[mask]) > 0:
                max_index = function(data[mask][:, 1])
                self.maxima = data[mask][max_index]
                self.plots["show maxima"].set_offsets(self.maxima)

        if len(self.clicks) == 2:
            self.max_index = function(data[mask][:, 1])
            plt.close(self.fig)

    def return_function(self):
        return np.array(self.maxima), self.max_index

    def run_condition(self):
        return self.plot_open


# %%
class LineFit(RangeSelectBase):
    def __init__(self, x, y, maxpoint=None):
        super().__init__(x, y)
        (self.plots["fit line"],) = self.ax.plot([], [], color="red")
        (self.plots["intercept"],) = self.ax.plot([], [], color="green")
        self.plots["text"] = self.ax.text(0, 0, "")

        self.maxpoint = maxpoint

        self.intercept_base = None
        self.params = None

    def custom_animation(self):
        data = self.data

        self._mouse_pointer()
        self._line_follow()
        self._range_highlight()

        mask = self.range_mask  # must be assigned with self._range_highlight()

        if (len(self.clicks) == 1) and (len(data[mask]) > 1):
            self.params, _ = curve_fit(line, data[mask][:, 0], data[mask][:, 1])
            print(r_squared(data[mask][:,1], data[mask][:,0], self.params[0], self.params[1]))
            
            p_x = line(data[mask][:,0], self.params[0], self.params[1])
            
            print(r2_score(data[mask][:,1], p_x))
            
            fit_x1 = data[mask][0][0]

            if self.maxpoint is None:
                fit_x2 = self.axes["xmax"]
            else:
                fit_x2 = self.maxpoint[0]

            fit_y1 = fit_x1 * self.params[0] + self.params[1]
            fit_y2 = fit_x2 * self.params[0] + self.params[1]

            self.plots["fit line"].set_data([fit_x1, fit_x2], [fit_y1, fit_y2])
            

            if self.maxpoint is not None:
                self.plots["intercept"].set_data(
                    [fit_x2, fit_x2], [fit_y2, self.maxpoint[1]]
                )
                self.plots["text"].set_text(str(self.maxpoint[1] - fit_y2))
                self.plots["text"].set_position([fit_x2, fit_y2])
                self.intercept_base = fit_y2

        if len(self.clicks) == 2:
            plt.close(self.fig)

    def return_function(self):
        return self.range_mask, self.intercept_base, self.params

    def run_condition(self):
        return self.plot_open


def line(x, a, b):
    return a * x + b

def r_squared(ydata, xdata, a, b):
    fitted = line(xdata, a, b)
    residuals = ydata - fitted
    #residuals = fitted - np.mean(ydata)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    #print(ss_res, ss_tot, ss_res / ss_tot)
    return 1-(ss_res / ss_tot)

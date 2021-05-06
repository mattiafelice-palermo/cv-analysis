#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 09:42:10 2021

@author: mattia
"""

import numpy as np
from scipy.optimize import curve_fit

import matplotlib

matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation

from PyQt5 import QtTest, QtCore


class CustomAnimation(animation.FuncAnimation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._running = False
        self._last_valute = None

    def new_frame_seq(self):
        # print("new_frame_seq")
        def internal_iterator():
            for e in animation.FuncAnimation.new_frame_seq(self):
                self._last_value = e
                if not self._running:
                    return self._last_value
                yield e

        return internal_iterator()

    def start(self):
        self._running = True
        self.event_source.start()

    def stop(self):
        self._running = False
        self.event_source.stop()


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        # Plot axes
        self.axes = fig.add_subplot(111)

        super(MplCanvas, self).__init__(fig)
        self.setStyleSheet("background-color:transparent;")
        self.figure.set_facecolor("none")
        self.figure.tight_layout(pad=0.1)

        # Reference to parent, if needed
        self.parent = None

        # For animations
        self.data = None
        self.bg_data = None
        self._ani = None
        self.animation = False
        self.plots = None
        self.return_data = None
        self.range_mask = None
        self.mouse_pos = [[0, 0], [0, 0]]
        self.clicks = []

        # Array containing True/False mask for the selected range
        self.range_mask = None

        # hide mouse cursor, needs a better implementation...
        # self.setCursor(QtCore.Qt.BlankCursor)

    def initialize(self, parent):
        self.parent = parent

    def cv_plot(self, filepath, cycle, checkboxes):
        # print(a_curr, a_ip, a_cf, c_curr, c_ip, c_cf)
        analysis = self.parent.file_handler.analyses[filepath][cycle]

        if self.animation is True:
            # dir(self._ani.event_source)
            # self._ani.event_source.stop()
            self.animation = False

        self.axes.cla()
        self.axes.set_xlabel("V vs ref [V]")
        self.axes.set_ylabel("i [mA]")

        if checkboxes.anode_curr:
            self.axes.plot(analysis.ox_voltage, analysis.ox_current, color="#1f77b4")
        if checkboxes.cathode_curr:
            self.axes.plot(analysis.red_voltage, analysis.red_current, color="#1f77b4")
        if checkboxes.anode_ip:
            pass
        if checkboxes.cathode_ip:
            pass
        if checkboxes.anode_capfit:
            pass
        if checkboxes.cathode_capfit:
            pass
        # self.figure.tight_layout(pad=0.1)
        self.draw()

    def animated_plot(self, data, bg_data, plot=None, other=None):
        self.data = np.vstack(data).T
        self.bg_data = np.vstack(bg_data).T
        self.plt_range = self._get_axes_range()  # get range from bg_x, bg_y

        # Clean figure and set new xlim, ylim
        # self.axes.cla()
        self.axes.set_xlim((self.plt_range["xmin"], self.plt_range["xmax"]))
        self.axes.set_ylim((self.plt_range["ymin"], self.plt_range["ymax"]))

        # runs animations
        self.plots = {
            "background": self.axes.plot([], [], color="#1f77b4")[0],
            "base": self.axes.plot([], [], color="#1f77b4")[0],
            "mouse": self.axes.scatter([], [], marker="+"),
            "line follow": self.axes.scatter([], []),
            "range hlight": self.axes.scatter([], []),
        }

        # Mode should indicate animated function to call in animation.FuncAnimation
        # create several self._animate variants to be called

        if plot == "maximum":
            function = self._maximum
            init_func = self._init_frame
        if plot == "minimum":
            function = self._minimum
            init_func = self._init_frame
        if plot == "linefit":
            function = lambda peak: self._linefit(0, other)
            init_func = self._init_frame
        if plot == "test":
            self.plots = {
                "test": self.axes.plot([], [])[0],
            }
            init_func = self._init_test
            function = self._animate_test

        # State that the animation is running
        self._ani = CustomAnimation(
            self.figure, function, init_func=init_func, interval=40, blit=True,
        )

        self._ani.start()

        # Update mouse position
        cid1 = self.figure.canvas.mpl_connect(
            "motion_notify_event", self._on_mouse_move
        )
        cid2 = self.figure.canvas.mpl_connect(
            "button_press_event", self._on_mouse_click
        )
        cid3 = self.figure.canvas.mpl_connect("resize_event", self._on_resize)

        self.loop = QtCore.QEventLoop()
        self.loop.exec()

        # self._ani.event_source.stop()
        self._ani.stop()

        # del self._ani
        # del self.plots
        # del self.loop

        self.animation = False

        return self.return_data

    def _init_test(self):
        return (self.plots["test"],)

    def _animate_test(self, i):
        print(f"Running animation with frame {i}")
        if len(self.clicks) == 1:
            self.clicks = []
            print(f"{i} before loop.quit()")
            self.loop.quit()
        print(f"{i} after loop.quit()")

        return tuple(self.plots.values())

    def _get_axes_range(self, c=0.05):
        xmin, xmax = np.amin(self.bg_data[:, 0]), np.amax(self.bg_data[:, 0])
        ymin, ymax = np.amin(self.bg_data[:, 1]), np.amax(self.bg_data[:, 1])

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
        self.plots["background"].set_data(self.bg_data[:, 0], self.bg_data[:, 1])
        return (self.plots["background"], self.plots["base"])

    def _maximum(self, i):
        print("in loop cause of _maximum")
        self.plots["show maximum"] = self.axes.scatter([], [], color="red")

        self._mouse_pointer()
        self._line_follow()
        self._range_highlight()

        mask = self.range_mask
        data = self.data

        if len(self.clicks) == 1:
            if len(data[mask]) > 0:
                max_index = np.argmax(data[mask][:, 1])
                maximum = data[mask][max_index]
                self.plots["show maximum"].set_offsets(maximum)

        if len(self.clicks) == 2:
            max_index = np.argmax(data[mask][:, 1])
            maximum = data[mask][max_index]
            self.return_data = maximum, max_index
            self.clicks = []
            print(f"{i} before loop.quit()")
            self.loop.quit()

            # should do all the savings here!

        print(f"{i} after loop.quit()")
        # return (self.plots["mouse"], self.plots["base"])
        return tuple(self.plots.values())

    def _minimum(self, i):
        print("in loop cause of _minimum")
        self.plots["show minimum"] = self.axes.scatter([], [], color="red")

        self._mouse_pointer()
        self._line_follow()
        self._range_highlight()

        mask = self.range_mask
        data = self.data

        if len(self.clicks) == 1:
            if len(data[mask]) > 0:
                min_index = np.argmin(data[mask][:, 1])
                minimum = data[mask][min_index]
                self.plots["show minimum"].set_offsets(minimum)

        if len(self.clicks) == 2:
            # self._ani.event_source.stop()

            min_index = np.argmin(data[mask][:, 1])
            minimum = data[mask][min_index]
            self.return_data = minimum, min_index
            self.clicks = []
            self.loop.quit()

            # should do all the savings here!

        # return (self.plots["mouse"], self.plots["base"])
        return tuple(self.plots.values())

    def _linefit(self, i, peak):
        print("in loop cause of _linefit")
        data = self.data

        (self.plots["fit line"],) = self.axes.plot([], [], color="red")
        (self.plots["intercept"],) = self.axes.plot([], [], color="green")
        self.plots["text"] = self.axes.text(0, 0, "")

        self._mouse_pointer()
        self._line_follow()
        self._range_highlight()

        mask = self.range_mask  # must be assigned with self._range_highlight()

        if (len(self.clicks) == 1) and (len(data[mask]) > 1):
            params, _ = curve_fit(line, data[mask][:, 0], data[mask][:, 1])
            # print(
            #    r_squared(
            #        data[mask][:, 1], data[mask][:, 0], self.params[0], self.params[1]
            #    )
            # )

            # p_x = line(data[mask][:,0], self.params[0], self.params[1])

            # print(r2_score(data[mask][:,1], p_x))

            fit_x1 = data[mask][0][0]

            if peak is None:
                fit_x2 = self.axes["xmax"]
            else:
                fit_x2 = peak[0]

            fit_y1 = fit_x1 * params[0] + params[1]
            fit_y2 = fit_x2 * params[0] + params[1]

            self.plots["fit line"].set_data([fit_x1, fit_x2], [fit_y1, fit_y2])

            if peak is not None:
                self.plots["intercept"].set_data([fit_x2, fit_x2], [fit_y2, peak[1]])
                self.plots["text"].set_text(str(peak[1] - fit_y2))
                self.plots["text"].set_position([fit_x2, fit_y2])
                intercept_base = fit_y2

        if len(self.clicks) == 2:
            # print(dir(self._ani.event_source))
            # print(dir(self._ani.event_source))
            self.clicks = []
            self.loop.quit()

        # return (self.plots["mouse"], self.plots["base"])
        return tuple(self.plots.values())

    def _on_mouse_move(self, event):
        # element 0 is the most recent mouse pos, element 1 is the previous
        if event.xdata is None:
            self.mouse_pos[0] = self.mouse_pos[1]
        else:
            self.mouse_pos[1] = self.mouse_pos[0]
            self.mouse_pos[0] = [event.xdata, event.ydata]

    def _on_mouse_click(self, event):
        # usability while zooming - to be improved with ESC catch?
        # if plt.get_current_fig_manager().toolbar.mode != "":
        #    return
        self.clicks.append([event.xdata, event.ydata])

    # def resizeEvent(self, event):
    #     if self.animation:
    #         self._ani.event_source.stop()
    #     else:
    #         super(MplCanvas, self).resizeEvent(event)

    def _on_resize(self, event):
        print("in _on_resize")
        print(self.animation)

        if self.animation:
            print("here!!")
            self._ani.event_source.stop()

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

    def _mouse_pointer(self):
        mouse = self.mouse_pos[0][0], self.mouse_pos[0][1]
        self.plots["mouse"].set_offsets([mouse[0], mouse[1]])


def line(x, a, b):
    return a * x + b

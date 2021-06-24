#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:25:35 2021

@author: mattia
"""
# from read_input import CyclicVoltammetry
# from cvanalysis import SCVAnalysis, CCVAnalysis
# from graphical_tools import MaxMinFinder, LineFit

import sys
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QAction, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5.Qt import QStandardItemModel, QStandardItem, QLayout

from MainWindow import Ui_MainWindow
from file_handler import FileHandler
from gui.GUI_dataclasses import MplChekboxStates


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.usrtype = (
            QtWidgets.QTreeWidgetItem.UserType
        )  # communication channel for data

        # Menubar
        menu = self.menubar

        file_menu = menu.addMenu("&File")
        add_folder_action = QAction(QIcon("icons/folder-plus.svg"), "Add folder", self)
        add_folder_action.triggered.connect(self.add_folder)
        file_menu.addAction(add_folder_action)

        # Manage files
        self.file_handler = FileHandler()

        # Selected CV
        self.selected_cv = None  # string containing CV identifier (filepath)

        # Add QTreeWidget and give it reference to parent to access file_handler
        # and MplCanvas
        self.CVTreeWidget.initialize(self)  # link parent to widget and set signals
        self.MplCanvas.initialize(self)
        # self.CVTreeWidget.show() # not useful

        # MplCanvas checkboxes
        self.anode_current_checkbox.stateChanged.connect(
            self.MplCanvas_checkbox_changed
        )
        self.anode_capfit_checkbox.stateChanged.connect(self.MplCanvas_checkbox_changed)
        self.anode_ip_checkbox.stateChanged.connect(self.MplCanvas_checkbox_changed)
        self.cathode_current_checkbox.stateChanged.connect(
            self.MplCanvas_checkbox_changed
        )
        self.cathode_capfit_checkbox.stateChanged.connect(
            self.MplCanvas_checkbox_changed
        )
        self.cathode_ip_checkbox.stateChanged.connect(self.MplCanvas_checkbox_changed)

        # Fit button
        self.manual_fit_button.clicked.connect(self.manual_fit)

    def add_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select folder")
        self.file_handler.add_folder(folder)
        self.CVTreeWidget.filetree_update()

    def on_CVTreeWidget_select(self, obj, col):
        parent = obj.parent()

        # Show SCVAnalysis data
        if parent is not None:
            cycle = obj.data(0, self.usrtype)
            cv = parent.data(0, self.usrtype)
            filepath = cv.filepath
            self.selected_cv = (filepath, cycle)

            analysis = self.file_handler.analyses[filepath][cycle]

            self.update_ipvalues()

            self.MplCanvas.cv_plot(filepath, cycle, self.MplCanvas_checkbox_states())

    def update_ipvalues(self):
        (filepath, cycle) = self.selected_cv
        analysis = self.file_handler.analyses[filepath][cycle]

        # Values of ip and voltage shown in main interface
        if analysis.anode_data.error is None:
            ipa_curr = analysis.anode_data.ip
            ipa_volt = analysis.anode_data.peak_volt
        else:
            ipa_curr, ipa_volt = "Automatic failed", "Automatic failed"

        if analysis.cathode_data.error is None:
            ipc_curr = analysis.cathode_data.ip
            ipc_volt = analysis.cathode_data.peak_volt
        else:
            ipc_curr, ipc_volt = "Automatic failed", "Automatic failed"

        self.anode_ip_value.setText(str(ipa_curr))
        self.anode_ippos_value.setText(str(ipa_volt))
        self.cathode_ip_value.setText(str(ipc_curr))
        self.cathode_ippos_value.setText(str(ipc_volt))

    def MplCanvas_checkbox_changed(self, s=None):
        if self.selected_cv is not None:
            checkboxes = self.MplCanvas_checkbox_states()

            (filepath, cycle) = self.selected_cv

            self.MplCanvas.cv_plot(filepath, cycle, checkboxes)

    def MplCanvas_checkbox_states(self):

        return MplChekboxStates(
            self.anode_current_checkbox.isChecked(),
            self.anode_ip_checkbox.isChecked(),
            self.anode_capfit_checkbox.isChecked(),
            self.cathode_current_checkbox.isChecked(),
            self.cathode_ip_checkbox.isChecked(),
            self.cathode_capfit_checkbox.isChecked(),
        )

    def manual_fit(self):
        # self.layout().setSizeConstraint(QLayout.SetFixedSize)

        checkboxes = self.MplCanvas_checkbox_states()
        (filepath, cycle) = self.selected_cv

        analysis = self.file_handler.analyses[filepath][cycle]

        if checkboxes.anode_curr and checkboxes.cathode_curr:
            bg_data = (analysis.voltage, analysis.current)
        elif checkboxes.anode_curr:
            bg_data = (analysis.ox_voltage, analysis.ox_current)
        elif checkboxes.cathode_curr:
            bg_data = (analysis.red_voltage, analysis.red_current)

        if checkboxes.anode_curr:
            max_data = self.MplCanvas.animated_plot(
                (analysis.ox_voltage, analysis.ox_current), bg_data, plot="maximum",
            )

            fit_data = self.MplCanvas.animated_plot(
                (analysis.ox_voltage, analysis.ox_current),
                bg_data,
                plot="linefit",
                other=max_data[0],
            )

            analysis.anode_data.ip = fit_data[1]
            analysis.anode_data.fit_mode = "manual"
            analysis.anode_data.error = None
            analysis.anode_data.peak_base = fit_data[0]
            analysis.anode_data.current_max = max_data[0][1]
            analysis.anode_data.peak_volt = max_data[0][0]
            analysis.anode_data.peak_index = max_data[1]
            analysis.anode_data.capacitive_fit = None
            analysis.anode_data.fit_data_bool = None  # to be implemented

            print(analysis.anode_data)

        self.update_ipvalues()

        if checkboxes.cathode_curr:
            max_data = self.MplCanvas.animated_plot(
                (analysis.red_voltage, analysis.red_current), bg_data, plot="minimum",
            )

            fit_data = self.MplCanvas.animated_plot(
                (analysis.red_voltage, analysis.red_current),
                bg_data,
                plot="linefit",
                other=max_data[0],
            )

            analysis.cathode_data.ip = fit_data[1]
            analysis.cathode_data.fit_mode = "manual"
            analysis.cathode_data.error = None
            analysis.cathode_data.peak_base = fit_data[0]
            analysis.cathode_data.current_max = max_data[0][1]
            analysis.cathode_data.peak_volt = max_data[0][0]
            analysis.cathode_data.peak_index = max_data[1]
            analysis.cathode_data.capacitive_fit = None
            analysis.cathode_data.fit_data_bool = None  # to be implemented

            print(analysis.cathode_data)

        self.update_ipvalues()

        self.MplCanvas.cv_plot(filepath, cycle, self.MplCanvas_checkbox_states())


app = QtWidgets.QApplication(sys.argv)

window = MainWindow()
window.show()

sys.exit(app.exec())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 12:30:19 2021

@author: mattia
"""

from PyQt5.QtWidgets import QTreeWidget, QTreeWidgetItem
import sys, os


class CVTreeWidget(QTreeWidget):
    def __init__(self, *args, **kwargs):
        super(QTreeWidget, self).__init__(*args, **kwargs)
        # print('I\'m widget!')
        self.parent = None  # Mainwindow
        self.file_handler = None  # for convenience? self.parent.file_handler...
        self.MplCanvas = None  # reference to plot widget

        self.usrtype = QTreeWidgetItem.UserType  # communication channel for data
        # self.itemClicked.connect(self.on_click)  # signal for clicking!

    def initialize(self, parent):
        self.parent = parent
        self.file_handler = parent.file_handler
        self.MplCanvas = parent.MplCanvas

        # Set signals!
        self.itemClicked.connect(self.parent.on_CVTreeWidget_select)

    def filetree_update(self):
        self.clear() # eliminate all items to avoid duplicating entries...
        items = []

        # index for QTreeWidgetItem, this is still bogus for me...
        usrtype = self.usrtype

        for filepath, cv in self.file_handler:
            filename = os.path.basename(cv.filepath)
            n_cycles = cv.settings["n_cycles"]
            file_item = QTreeWidgetItem([filename])
            file_item.setData(0, usrtype, cv)  # usrtype sigh...

            for i in range(n_cycles):
                child = QTreeWidgetItem([str(i)])
                child.setData(0, usrtype, i)
                file_item.addChild(child)

            items.append(file_item)

        self.insertTopLevelItems(0, items)

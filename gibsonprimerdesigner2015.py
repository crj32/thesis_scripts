#!/usr/bin/python

# Gibson Vector Designer 2015, an open source project for biologists

import sys
import re
import time
import os
from Bio import SeqIO
import sqlite3 as lite
from PySide.QtCore import *
from PySide import QtGui
from PySide import QtCore
import collections
import itertools

import inputmessage

# Function for doing reverse complement

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases 

# Globals

gibson_dict = {}
gibson_sequence_names = []
global_parts = collections.OrderedDict()

# Class for the graphical user interface
        
class Ui_Gibson(object):
    def setupUi(self, Gibson):
        Gibson.setObjectName("Gibson")
        Gibson.resize(1061, 667)
        self.centralwidget = QtGui.QWidget(Gibson)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(40, 60, 361, 66))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_2 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.label = QtGui.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.lineEdit_2 = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.gridLayout.addWidget(self.lineEdit_2, 1, 1, 1, 1)
        self.lineEdit = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout.addWidget(self.lineEdit, 0, 1, 1, 1)
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(180, 20, 66, 21))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(520, 20, 121, 21))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtGui.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(150, 200, 131, 21))
        self.label_5.setObjectName("label_5")
        self.listWidget = QtGui.QListWidget(self.centralwidget)
        self.listWidget.setGeometry(QtCore.QRect(450, 61, 251, 421))
        self.listWidget.setObjectName("listWidget")
        self.gridLayoutWidget_2 = QtGui.QWidget(self.centralwidget)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(40, 240, 361, 151))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_2 = QtGui.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_7 = QtGui.QLabel(self.gridLayoutWidget_2)
        self.label_7.setObjectName("label_7")
        self.gridLayout_2.addWidget(self.label_7, 1, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem, 4, 0, 1, 1)
        self.pushButton_2 = QtGui.QPushButton(self.gridLayoutWidget_2)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout_2.addWidget(self.pushButton_2, 5, 0, 1, 1)
        self.pushButton_10 = QtGui.QPushButton(self.gridLayoutWidget_2)
        self.pushButton_10.setObjectName("pushButton_10")
        self.gridLayout_2.addWidget(self.pushButton_10, 5, 1, 1, 1)
        self.label_8 = QtGui.QLabel(self.gridLayoutWidget_2)
        self.label_8.setObjectName("label_8")
        self.gridLayout_2.addWidget(self.label_8, 2, 0, 1, 1)
        self.label_9 = QtGui.QLabel(self.gridLayoutWidget_2)
        self.label_9.setObjectName("label_9")
        self.gridLayout_2.addWidget(self.label_9, 3, 0, 1, 1)
        self.lineEdit_4 = QtGui.QLineEdit(self.gridLayoutWidget_2)
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.gridLayout_2.addWidget(self.lineEdit_4, 2, 1, 1, 1)
        self.lineEdit_3 = QtGui.QLineEdit(self.gridLayoutWidget_2)
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.gridLayout_2.addWidget(self.lineEdit_3, 3, 1, 1, 1)
        self.lineEdit_5 = QtGui.QLineEdit(self.gridLayoutWidget_2)
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.gridLayout_2.addWidget(self.lineEdit_5, 1, 1, 1, 1)
        self.horizontalLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(40, 130, 361, 51))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pushButton = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout.addWidget(self.pushButton)
        self.pushButton_6 = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.pushButton_6.setObjectName("pushButton_6")
        self.horizontalLayout.addWidget(self.pushButton_6)
        self.listWidget_2 = QtGui.QListWidget(self.centralwidget)
        self.listWidget_2.setGeometry(QtCore.QRect(760, 60, 251, 421))
        self.listWidget_2.setObjectName("listWidget_2")
        self.label_6 = QtGui.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(830, 10, 111, 41))
        self.label_6.setObjectName("label_6")
        self.verticalLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(760, 500, 251, 95))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.pushButton_3 = QtGui.QPushButton(self.verticalLayoutWidget)
        self.pushButton_3.setObjectName("pushButton_3")
        self.verticalLayout.addWidget(self.pushButton_3)
        self.pushButton_4 = QtGui.QPushButton(self.verticalLayoutWidget)
        self.pushButton_4.setObjectName("pushButton_4")
        self.verticalLayout.addWidget(self.pushButton_4)
        self.pushButton_9 = QtGui.QPushButton(self.verticalLayoutWidget)
        self.pushButton_9.setObjectName("pushButton_9")
        self.verticalLayout.addWidget(self.pushButton_9)
        self.verticalLayoutWidget_2 = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(450, 500, 251, 95))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.pushButton_8 = QtGui.QPushButton(self.verticalLayoutWidget_2)
        self.pushButton_8.setObjectName("pushButton_8")
        self.verticalLayout_2.addWidget(self.pushButton_8)
        self.pushButton_7 = QtGui.QPushButton(self.verticalLayoutWidget_2)
        self.pushButton_7.setObjectName("pushButton_7")
        self.verticalLayout_2.addWidget(self.pushButton_7)
        self.pushButton_5 = QtGui.QPushButton(self.verticalLayoutWidget_2)
        self.pushButton_5.setObjectName("pushButton_5")
        self.verticalLayout_2.addWidget(self.pushButton_5)
        self.textBrowser = QtGui.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(40, 440, 361, 151))
        self.textBrowser.setObjectName("textBrowser")
        self.label_11 = QtGui.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(170, 410, 81, 21))
        self.label_11.setObjectName("label_11")
        Gibson.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(Gibson)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1061, 23))
        self.menubar.setObjectName("menubar")
        Gibson.setMenuBar(self.menubar)

        QtCore.QMetaObject.connectSlotsByName(Gibson)

        Gibson.setWindowTitle(QtGui.QApplication.translate("Gibson", "Gibson Primer Designer", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Gibson", "Part Name:", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Gibson", "Sequence: ", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Gibson", "<html><head/><body><p><img src=\":/newPrefix/Inputs.png\"/></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("Gibson", "<html><head/><body><p><img src=\":/newPrefix/Inventory.png\"/></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("Gibson", "<html><head/><body><p><img src=\":/newPrefix/parameters.png\"/></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("Gibson", "<html><head/><body><p>Tm (<span style=\" vertical-align:super;\">o</span>C)</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_2.setText(QtGui.QApplication.translate("Gibson", "Build Vector Primers", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_10.setText(QtGui.QApplication.translate("Gibson", "View vector sequence", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("Gibson", "Min primer  length (bp)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("Gibson", "Max primer length (bp)", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Gibson", "Add Part", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_6.setText(QtGui.QApplication.translate("Gibson", "Load from file", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("Gibson", "<html><head/><body><p><img src=\":/newPrefix/database.png\"/></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_3.setText(QtGui.QApplication.translate("Gibson", "Transfer part to inventory", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_4.setText(QtGui.QApplication.translate("Gibson", "View sequence of part", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_9.setText(QtGui.QApplication.translate("Gibson", "Delete part from database", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_8.setText(QtGui.QApplication.translate("Gibson", "Transfer part to database", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_7.setText(QtGui.QApplication.translate("Gibson", "View sequence of part", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_5.setText(QtGui.QApplication.translate("Gibson", "Reset vector parts", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("Gibson", "<html><head/><body><p><img src=\":/newPrefix/System.png\"/></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        
        # Initialise objects to hold data

        self.number_of_parts = None
        self.size_of_vector = None

        # These are two lists for the primer building algorithm

        self.sequences = []
        self.sequences_names = []

        # This is a dictionary for the inventory name&sequence to hold them for transfer to DB

        self.inventory = collections.OrderedDict()

        # Connect everything up

        self.pushButton.clicked.connect(self.addpart)
        self.pushButton_5.clicked.connect(self.clearparts)
        self.pushButton_6.clicked.connect(self.loadpart)
        self.pushButton_2.clicked.connect(self.buildprimers)
        self.pushButton_8.clicked.connect(self.addtodict)
        self.pushButton_3.clicked.connect(self.exportfromdb)
        self.pushButton_9.clicked.connect(self.deletefromdb)
        self.pushButton_10.clicked.connect(self.view_vector_sequence)

        # For toolbar

        m = 45

        # Initialise the database

        con = lite.connect('./test.db')

        l = []

        with con:    
            cur = con.cursor()    
            cur.execute("SELECT * FROM Seqs")
            rows = cur.fetchall()
            for row in rows:
                for item in row:
                    l.append(item)
        
        print l

        for i in range(0,len(l),2):
            x = l[i]
            self.listWidget_2.addItem(x)

        ######################## For functionality (non database)

    def quit(self):
        sys.exit()

    def addpart(self):

        toaddname = self.lineEdit.text()
        toaddsequence = self.lineEdit_2.text()

        if str(toaddsequence) == "":
            self.textBrowser.append("error, you cant join an empty part")
        else:
            # to add to the list
            self.sequences.append(str(toaddsequence)) 
            self.sequences_names.append(str(toaddname)) # add the names back in later on for primer visualisation
            self.listWidget.addItems([toaddname])
            self.textBrowser.append("parts added")
            self.update_statistics()
            # to add to the dictionary
            self.inventory[str(toaddname)] = str(toaddsequence)
            print self.inventory

    def clearparts(self):

        # for list
        self.listWidget.clear()
        self.sequences = []
        self.textBrowser.append("parts cleared")
        self.update_statistics()
        # for dict
        self.inventory.clear()

    def loadpart(self): # test this, also upgrade so user can input multi .fasta file

        path, _ = QtGui.QFileDialog.getOpenFileName(mySW, "Open File", os.getcwd())
        
        inputfasta = SeqIO.parse(path, 'fasta')
        for record in inputfasta:
            # for list
            self.sequences.append(record.seq.tostring())
            self.sequences_names.append(record.id)
            # for dict
            self.inventory[record.id] = record.seq.tostring()

            self.listWidget.addItems([record.id])

        self.update_statistics()

        print self.sequences   
        print self.sequences_names 

    def update_statistics(self):

        self.number_of_parts = self.listWidget.count()
        self.textBrowser.append("no. of parts: " + str(self.number_of_parts))

        vector = ''.join(self.sequences)

        self.size_of_vector = str(len(vector)/1000) + " Kbp"
        self.textBrowser.append("vector size: " + self.size_of_vector)

    def view_vector_sequence(self):

        global global_parts
        global_parts = self.inventory

        self.vector_window = vector_window()
        self.vector.window.show()

    ######################## Database functions

    def addtodict(self):

        userselection1 = self.listWidget.currentItem()
        userselection = userselection1.text()
        print userselection
        seq4db =  {}
        for k, v in self.inventory.iteritems():
            if k == userselection: # replace with user selection
                seq4db[k] = v
        seq4db2 = [[k, v] for k, v in seq4db.iteritems()]

        con = lite.connect('./test.db')

        # Add to database

        with con:
            cur = con.cursor()    
            #cur.execute("DROP TABLE IF EXISTS Seqs")
            cur.execute("CREATE TABLE IF NOT EXISTS Seqs(Id TEXT, Seq TEXT)")
            cur.executemany("INSERT INTO Seqs VALUES(?, ?)", seq4db2)

        # Display in list in GUI

        for k in seq4db.keys():
            self.listWidget_2.addItems([k])

    def exportfromdb(self):

        con = lite.connect('./test.db')

        userselection1 = self.listWidget_2.currentItem()
        userselection = userselection1.text()
    
        var = [userselection] # replace with user selected input from database list
        l = []
        for item in var:
            with con:    
                cur = con.cursor()    
                cur.execute("SELECT * FROM Seqs WHERE Id=?", (item,))
                rows = cur.fetchall()
                for row in rows:
                    for item in row:
                        l.append(item)               

        print self.inventory

        self.inventory[l[0].encode("utf-8")] = l[1].encode("utf-8")

        self.listWidget.addItem(l[0])

        print self.inventory

    def deletefromdb(self):

        con = lite.connect('./test.db')

        userselection1 = self.listWidget_2.currentItem()
        userselection = userselection1.text()

        with con:
            cur = con.cursor()
            cur.execute("DELETE FROM Seqs WHERE Id=?", (userselection,))

        listItems = self.listWidget_2.selectedItems()
        for item in listItems:
            self.listWidget_2.takeItem(self.listWidget_2.row(item))

    ###################### Gibson primer selection algorithm 

    def buildprimers(self):

        # Re-define self.sequences and self.sequences_names from self.inventory
        # You can remove all non dictionary code earlier on at somepoint

        print self.sequences_names
        print self.sequences

        self.sequences_names = []
        self.sequences = []

        for k in self.inventory.keys():
            #print k
            self.sequences_names.append(k)

        for k, v in self.inventory.items():
            #print k, v
            self.sequences.append(v)

        #

        ERROR = 0 # set the error code to control the output window

        self.textBrowser.append("building gibson primers")

        # define possible errors

        if (self.lineEdit_5.text() == "") or (self.lineEdit_4.text() == "") or (self.lineEdit_3.text()) == "":
            ERROR = 1
            self.textBrowser.append("error, parameters were not set")

        self.tm = int(self.lineEdit_5.text())
        self.min_size = int(self.lineEdit_4.text())
        self.max_size = int(self.lineEdit_3.text())

        # Creating two lists of sequence IDs

        seqnumber = len(self.sequences)

        self.sequences1 = []
        self.sequences2 = []

        for i in range(1,seqnumber+1,1):
            theid = "sequence" + str(i)
            self.sequences1.append(theid)
            self.sequences2.append(theid)

        self.sequences1 = self.sequences1[:seqnumber-1]
        self.sequences2 = self.sequences2[1:seqnumber+1]

        #print self.sequences1
        #print self.sequences2

        dict = collections.OrderedDict() # this creates a dictionary that maintains the initial order
        m = 1

        # Get primer sequences 5' and 3' without overlaps

        for sequence in self.sequences:
            n = self.min_size
            n2 = self.min_size
            primer = sequence[:n]
            primer2 = sequence[-n:]
            Tm = 0
            Tm2 = 0
            while (Tm < self.tm):
                n += 1
                primer = sequence[:n]
                Tm = 64.9 + 41*(primer.count("G") + primer.count("C") - 16.4)/(primer.count("A") + primer.count("T") + primer.count("G") + primer.count("C"))
                if n >= self.max_size:
                    print "warning maximum size reached for forward primer part", m
                    break
            while (Tm2 < self.tm):
                n2 += 1
                primer2 = sequence[-n2:]
                Tm2 = 64.9 + 41*(primer2.count("G") + primer2.count("C") - 16.4)/(primer2.count("A") + primer2.count("T") + primer2.count("G") + primer2.count("C"))
                if n >= self.max_size:
                    print "warning maximum size reached for reverse primer part", m
                    break
            if (Tm and Tm2) >= self.tm and (n and n2 < self.max_size):
                x = "sequence" + str(m)
                m += 1
                fiveprime = sequence[:20]
                threeprime = sequence[-20:]
                dict[x] = [primer, primer2, Tm, Tm2, fiveprime, threeprime]
            else:
                self.textBrowser.append("cannot find primers within the parameter range")
                ERROR = 1
                
        # 5' ENDS extract and integrate

        ends = {}
        counter = 2

        for key, value in sorted(dict.items()):
            if key in self.sequences1:
                name = "sequence" + str(counter)
                ends[name] = [value[5]] 
                counter += 1
            else:
                name = "sequence" + str(1)
                ends[name] = [value[5]]

        for k, v in sorted(dict.items()):
            for k2, v2 in sorted(ends.items()):
                if k == k2:
                    dict[k] = [v2[0]+v[0], v[1], v[2], v[3], v[4], v[5]] # modify dict with the 5' ends

        # 3' ENDS extract and integrate

        ends2 = {}
        counter2 = 1

        for key, value in sorted(dict.items()): # This loop is not correct
            if key in self.sequences2:
            # for everything apart from sequence 1
                name = "sequence" + str(counter2)
                ends2[name] = [value[4]] 
                counter2 += 1
            else:
                name = "sequence" + str(seqnumber) # This variable will alter based on number of sequences added
                ends2[name] = [value[4]]

        for k, v in sorted(dict.items()):
            for k2, v2 in sorted(ends2.items()):
                if k == k2:
                    dict[k] = [v[0], v[1]+v2[0], v[2], v[3], v[4], v[5]]    

        # Add a loop that reverse complements the second primer

        for k, v in sorted(dict.items()):
            dict[k] = [v[0], reverse_complement(v[1]), v[2], v[3], v[4], v[5]] 

        # Put back in the original sequence names into the dict or newdict  

        i = 0

        for k in dict.keys():
            newid = self.sequences_names[i]
            dict[newid] = dict.pop(k)
            i += 1

        # Add in primer lengths

        for k, v in sorted(dict.items()):
            lenforward = len(v[0]) - 20
            lenreverse = len(v[1]) - 20
            dict[k] = [v[0], v[1], v[2], v[3], v[4], v[5], lenforward, lenreverse]

        # print out results in new window in .csv format with a header

        global gibson_dict
        gibson_dict = dict # I could not get the other class to recognise the dictionary any other way

        global gibson_sequence_names
        gibson_sequence_names = self.sequences_names

        if ERROR == 0:
            self.gibsonoutputwindow = gibsonoutputwindow()
            self.gibsonoutputwindow.show()
            self.textBrowser.append("primers built!")

class ControlGibsonWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):
        super(ControlGibsonWindow, self).__init__(parent)
        self.ui = Ui_Gibson()
        self.ui.setupUi(self)

        self.show()

##################################### OUTPUT WINDOWS (I dont think you can import these as modules?)

class gibsonoutputwindow(QtGui.QMainWindow): 

        def __init__(self, parent=None):
            super(gibsonoutputwindow, self).__init__(parent)
            self.setWindowTitle('Gibson primer builder output')

            self.setGeometry(500, 500, 1500, 500)
            self.grid = QtGui.QGridLayout()
            self.grid.setSpacing(10)

            self.gibson_output = QtGui.QTextEdit(self)
            self.grid.addWidget(self.gibson_output, 0, 0)

            self.toolbar = self.addToolBar("Quit")

            SaveAction = QtGui.QAction('Save as .csv file', self)
            SaveAction.triggered.connect(self.save)

            self.toolbar.addAction(SaveAction)

            self.mainWidget = QtGui.QWidget() 
            self.mainWidget.setLayout(self.grid)
            self.setCentralWidget(self.mainWidget)

            self.show()

            # change text to include gibson primer output # how to get dictionary from other class?

            global gibson_dict
            print gibson_dict
            global gibson_sequence_names
            print gibson_sequence_names

            gibson_list = []

            for k, v in gibson_dict.items():
                gibson_list.append(k)
                gibson_list.extend((v[0], v[1], str(v[2]), str(v[3]), str(v[6]), str(v[7])))

            print gibson_list

            empty = "part name, primer forward, primer reverse, tm forward, tm reverse, length forward, length reverse\n"
            for i in range(0, len(gibson_list),7):
                empty += ','.join(gibson_list[i:i+7]) + "\n"

            self.gibson_output.setText(empty)

            self.towrite = empty

        def save(self):

            path2, _ = QtGui.QFileDialog.getSaveFileName(self, "Save file", "", ".csv")
            towrite = self.towrite
            with open(path2+".csv", "w") as text_file:
                text_file.write(towrite)

class vector_window(QtGui.QMainWindow):

        def __init__(self, parent=None):
            super(vector_window, self).__init__(parent)
            self.setWindowTitle('Vector sequence viewer')

            self.setGeometry(500, 500, 500, 500)
            self.grid = QtGui.QGridLayout()
            self.grid.setSpacing(10)

            self.vector_viewer = QtGui.QTextEdit(self)
            self.grid.addWidget(self.vector_viewer, 0, 0)

            self.mainWidget = QtGui.QWidget() 
            self.mainWidget.setLayout(self.grid)
            self.setCentralWidget(self.mainWidget)

            self.show()

            # print out vector sequence

            list = []
            global global_parts
            for k,v in global_parts.iteritems():
                list.append(v)
            string = ''.join(list)
            self.vector_viewer.setText(string)    

############################################## End of classes

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    mySW = ControlGibsonWindow()
    mySW.show()
    sys.exit(app.exec_())

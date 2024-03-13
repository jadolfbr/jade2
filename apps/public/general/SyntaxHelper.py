#!/usr/bin/env python

import os
import re
from tk import *
from tkinter import *
from tkinter import filedialog as tkFileDialog
from tkinter import simpledialog as tkSimpleDialog
import pickle

class GUI():
	def __init__(self, main):
		'''
		Initializes Tk and all variables
		'''
		self.main = main
		self.pwd = (self.location())[0]
		self.SettingsPath = "SyntaxData.bin" #DEFAULT - If you load a text file - that text will convert to a bin, and then you can modify it from there.
		self.language = StringVar()
		self.language.set("")
		self.typeEntry = StringVar() #Type you enter in
		#Data Loading.
		self.typeSelect = StringVar() #Acts as a place holder for types.
		self.codeSelect = StringVar() #Acts as a place holder for code
		self.SynData = dict() #Main Data Structure [language][type][code]=description
		self.SynData = self.checkSettings()
		if len(self.SynData)==0:
			self.language.set("Languages")
		else:

			OPTIONS=list(self.SynData.keys())
			OPTIONS.sort()
			self.language.set(OPTIONS[0])

		#TK
		self.languagesTk = OptionMenu(self.main, self.language, self.language.get())
		self.languagesTk2 = OptionMenu(self.main, self.language, list(self.SynData.keys()))
		#Conditional to load out of the dictionary:
		if len(self.SynData)!=0:
			self.get_languages_from_dictionary_and_set()

		self.language.trace_variable('w', self.updateTypes)
		self.typesTk = Listbox()
		self.codeTk = Listbox()
		self.frameHelp = Frame(self.main, bd=3, relief=GROOVE);
		self.textHelp = Text(self.frameHelp,wrap="word", height = 20, width = 50, background = 'white')
		self.MenBar = Menu(self.main)
		self.MenInsert = Menu(self.main, tearoff=0);
		self.MenInsert.add_command(label="Save", command = lambda: self.writeData())
		self.MenInsert.add_command(label="New Language", command=lambda: self.importLanguage())
		#self.MenInsert.add_command(label = "Import Text File")
		self.MenInsert.add_command(label = "Import Binary", command = lambda: self.importBinary())
		self.MenBar.add_cascade(label = "Import/Insert", menu=self.MenInsert)

		self.MenExport = Menu(self.main, tearoff=0);
		#self.MenExport.add_command(label="Export Text File", command=lambda:self.exportText())
		self.MenExport.add_command(label="Export Binary File", command = lambda: self.exportBinary())
		self.MenBar.add_cascade(label = "Export", menu=self.MenExport)
		self.MenEdit = Menu(self.main, tearoff=0);
		self.MenEdit.add_command(label = "Modify Code Entry", command=lambda: self.modifyCode())
		self.MenEdit.add_command(label = "Modify Description Entry", command=lambda: self.modifyText())
		self.MenEdit.add_separator()
		self.MenEdit.add_command(label = "Delete Code Entry", command = lambda: self.deleteCodeEntry())
		self.MenEdit.add_command(label = "Delete Type", command = lambda: self.deleteType())
		self.MenEdit.add_command(label = "Delete Language", command = lambda: self.deleteLanguage())
		self.MenBar.add_cascade(label = "Edit", menu=self.MenEdit)
		#Tk Insert New stuff

		self.typeEntryTk = Entry(self.main, textvariable=self.typeEntry)
		self.typeLab = Label(self.main, text = "Type")
		self.frameCode = Frame(self.main, bd=3, relief=GROOVE); self.textCode = Text(self.frameCode, wrap="word", height=3, width=35, background='white')
		self.CodeLab = Label(self.main, text ="Code")
		self.frameDesc = Frame(self.main, bd=3, relief=GROOVE); self.textDesc = Text(self.frameDesc, wrap="word", height=5, width=30, background='white')
		self.DescLab = Label(self.main, text = "Description")

		self.Add = Button(self.main, text = "Add Syntax Note", command = lambda: self.addSyntax_and_update())

		#Bindings
		self.typesTk.bind("<ButtonRelease-1>", lambda event: self.updateCodes())
		self.codeTk.bind("<ButtonRelease-1>", lambda event:  self.insertText())

		if self.language.get()!="Languages":
			self.updateTypes(1, 2, 3)


	def gridTk(self, row, column):
		'''
		Grids Tk
		'''


		#Get info about what to do:

		self.codeTk.config(width = 40)

		self.languagesTk.grid(row = row, column=column, sticky=W+E)
		self.typesTk.grid(row=row+1, column = column+1, sticky=W+E)
		self.codeTk.grid(row=row+1, column=column+2, sticky=W+E)
		self.textHelp.grid(row=row, column=column+3, columnspan = 2, rowspan = 2, sticky = W+E)
		self.frameHelp.grid(row = row, column=column+3, columnspan=2, rowspan = 2, sticky = W+E, padx=6, pady=7)

		#Entry
		self.languagesTk2.grid(row=row+2, column=column)
		self.typeEntryTk.grid(row = row+2, column=column+1)
		self.frameCode.grid(column=column+2, row = row+2); self.textCode.grid(column=column+3, row=row+2);
		self.frameDesc.grid(column=column+3, row=row+2); self.textDesc.grid(column=column+4, row=row+2)

		self.typeLab.grid(column=column+1, row=row+3)
		self.CodeLab.grid(column=column+2, row=row+3)
		self.DescLab.grid(column=column+3, row=row+3)

		self.Add.grid(column=column+4, row = row+2)
		self.main.config(menu=self.MenBar)


####################Commands.  Tried for them not to get rediculous, but they did anyway.
	def checkSettings(self):
		'''
		Checks if data exists.  Loads the binary data.
		'''
		if os.path.exists(self.SettingsPath):
			SynData = pickle.load(open(self.SettingsPath))
		else:
			SynData = dict()
		return SynData

	def writeData(self):
		'''
		Dumps SynData as a binary using pickle.
		'''
		pickle.dump(self.SynData, open(self.SettingsPath, 'wb'))

	def importLanguage(self):
		'''
		Puts a new language into the dictionary.  Saves binary.
		'''

		lang = tkSimpleDialog.askstring(title = "Language", prompt = "Please Enter the new language")
		self.SynData[lang]=dict()
		self.get_languages_from_dictionary_and_set()
		#self.writeData()



	def get_languages_from_dictionary_and_set(self):
		'''
		Sets up the GUI to display info from loaded binary dictionary.
		'''

		OPTIONS = []
		print(self.SynData)
		for i in self.SynData: OPTIONS.append(i); OPTIONS.sort()
		print(OPTIONS)
		if OPTIONS.__contains__("Languages"):
			self.SynData.pop("Languages")
		self.updateList(self.languagesTk, OPTIONS)
		self.updateList(self.languagesTk2, OPTIONS)

	def updateList(self, optmenu, options):
		optmenu["menu"].delete(0, END)
		for i in options:
				optmenu["menu"].add_command(label=i, command=lambda temp = i: optmenu.setvar(optmenu.cget("textvariable"), value = temp))

	def addSyntax_and_update(self):
		'''
		What happens when you click 'add syntax'
		Adds the data into dictionary.  Saves binary.  Updates GUI.
		'''

		type = self.typeEntry.get()
		code = self.textCode.get(1.0, END); self.textCode.delete(1.0, END)
		text = self.textDesc.get(1.0, END); self.textDesc.delete(1.0, END)
		#First, we add to an existing type or create a new one.
		if type in self.SynData[self.language.get()]:
			self.SynData[self.language.get()][type][code]=text
		else:
			self.SynData[self.language.get()][type]=dict()
			self.SynData[self.language.get()][type][code]=text
		self.updateTypes(1,2,3)
		#self.writeData()

	def updateTypes(self, name, index, mode):
		'''
		Updates the type ListBox to reflect new info.
		'''

		lang = self.language.get()
		self.typesTk.delete(0, END)
		TYPES = list(self.SynData[lang].keys())
		TYPES.sort()
		self.typesTk.insert(END, *TYPES)
		self.codeTk.delete(0, END)
		self.textHelp.delete(1.0, END)

	def updateCodes(self, type=0):
		'''
		Updates the code ListBox to reflect new info.
		'''

		if type==0:
			CODES = list(self.SynData[self.language.get()][self.typesTk.get(self.typesTk.curselection())].keys())
			self.typeSelect.set(self.typesTk.get(self.typesTk.curselection()))
		else:
			CODES = list(self.SynData[self.language.get()][self.typeSelect.get()].keys())

		CODES.sort()
		self.codeTk.delete(0, END)
		self.codeTk.insert(END, *CODES)

		if type==0:
			self.typeEntry.set(self.typesTk.get(self.typesTk.curselection()))
		else:
			pass

	def insertText(self):
		'''
		Inserts the description into the TextBox frame.
		'''

		self.textHelp.delete(1.0, END)
		#self.textHelp.insert(1.0, self.codeTk.get(self.codeTk.curselection())+"\n\n")
		self.textHelp.insert(1.0, self.SynData[self.language.get()][self.typeSelect.get()][self.codeTk.get(self.codeTk.curselection())]+"\n\n")
		self.codeSelect.set(self.codeTk.get(self.codeTk.curselection()))
	def deleteCodeEntry(self):
		'''
		Deletes the selected code entry.
		'''

		self.SynData[self.language.get()][self.typeSelect.get()].pop(self.codeTk.get(self.codeTk.curselection()))
		self.updateCodes(type=1)
		#self.writeData()
	def deleteType(self):
		'''
		Deletes all codes of a specific type.
		'''

		self.SynData[self.language.get()].pop(self.typesTk.get(self.typesTk.curselection()))
		self.updateTypes(1,2,3)
		#self.writeData()
	def deleteLanguage(self):
		'''
		Deletes all data of the language.
		'''

		self.SynData.pop(self.language.get())
		self.get_languages_from_dictionary_and_set()
		#self.writeData()

	def modifyCode(self):
		'''
		Modifies the existing code through a dialog box.
		'''

		newCode = tkSimpleDialog.askstring(title="Modify Code", prompt="New Code", initialvalue=self.codeSelect.get())
		oldText = self.SynData[self.language.get()][self.typeSelect.get()][self.codeSelect.get()]
		self.SynData[self.language.get()][self.typeSelect.get()].pop(self.codeSelect.get())
		self.SynData[self.language.get()][self.typeSelect.get()][newCode]=oldText
		self.updateCodes(1)
		#self.writeData()

	def modifyText(self):
		'''
		Modifies the existing text by parsing the textbox and saving the modified description.
		'''

		newEntry = self.textHelp.get(1.0, END)
		self.SynData[self.language.get()][self.typeSelect.get()][self.codeSelect.get()]=newEntry
		print("Entry Updated.")
		#self.writeData()

	def exportText(self):
		'''
		Will be used to export all the data to a nicely formatted text file.
		'''

		pass
	def exportBinary(self):
		'''
		Exports the dictionary as a binary so you can share data you have accumulated.
		'''

		FILE = tkFileDialog.asksaveasfile(initialdir = self.pwd, mode='wb')
		pickle.dump(self.SynData, FILE)
	def importText(self):
		pass
	def importBinary(self):
		'''
		Imports the binary.  Replaces what you were working on with the new binary.
		'''

		self.textHelp.delete(1.0, END)
		self.SettingsPath = tkFileDialog.askopenfilename(initialdir = self.pwd)

		self.SynData.clear()
		self.SynData = pickle.load(open(self.SettingsPath))
		OPTIONS=list(self.SynData.keys())
		OPTIONS.sort()
		self.language.set(OPTIONS[0])
		if len(self.SynData)!=0:
			self.get_languages_from_dictionary_and_set()

		self.SettingsPath = self.pwd +"/temp.bin"
		return self.SynData

	def location(self):
		'''
		Allows the script to be self-aware of it's path.
		So that it can be imported from anywhere.
		'''

		p = os.path.abspath(__file__)
		pathSP = os.path.split(p)
		return pathSP

if __name__ == '__main__':
	'''
	If the program is called directly, by 'python SyntaxHelper.py', this runs.
	If it is imported into another program, this does not.
	'''

	mainWin = Tk()
	mainWin.title("Syntax Helper")
	GuiObjects = GUI(mainWin)
	GuiObjects.gridTk(0, 0)
	mainWin.mainloop()


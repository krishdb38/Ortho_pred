# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:21:35 2020

@author: user
"""


import tkinter as tk

def setfullname():
    f_name = first_entry.get()
    s_name = surname_entry.get()
    fullname.set(f_name+" "+s_name)

root = tk.Tk()
fullname = tk.StringVar()

fullname.set("test")

first_label = tk.Label(root, text = "first name")
first_entry = tk.Entry(root)

surnameLabel = tk.Label(root,text = "Surname")
surname_entry = tk.Entry(root,)

fullname_label = tk.Label(root, text = "Full Name")
fullname_entry = tk.Entry(root, textvariable = fullname)

btn_Fullname = tk.Button(root,text = "show_full_name",command = setfullname)




first_label.grid(row =0, column = 0, padx = 15, pady = 15)
first_entry.grid(row = 0, column =1, padx = 15, pady = 15)

surnameLabel.grid(row = 1, column =0, padx =15 , pady = 15)
surname_entry.grid(row =1 , column = 1, padx =15, pady = 15)

fullname_label.grid(row = 2, column = 0 , padx = 15, pady = 15)
fullname_entry.grid(row = 2, column =1, padx = 15, pady = 15)

btn_Fullname.grid(columnspan = 2, padx = 15 , pady = 15)


root.mainloop()
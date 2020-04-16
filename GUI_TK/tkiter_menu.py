# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:34:57 2020

Example showing for tkinter and ttk how to:
    --1.Make a menubar with menu's
    --2.Put menu items on the menu's
    3. Establish callback functions for the menu items, that is,
    functions that are called when a menu item is selected.
    

@author: Krishna
"""

import tkinter 
from tkinter import ttk

class Data():
    def __init__(self):
        self.number = 0
        self.number_label = None

def main():
    data = Data()
    root = tkinter.Tk()
    
    main_frame = ttk.Frame(root, padding = 50)
    main_frame.grid()
    
    label_text = "The number is " + str(data.number)
    label = ttk.Label(main_frame, text = label_text)
    label.grid()
    
    data.number_label = label
    
    # The default is for menus to be "tear-off" --they can be dragged
    # off the menubar. Use whichever style best your GUI
    
    root.option_add("*tearoff",False)
    
    # ! step 1 Make the menu bar
    menubar = tkinter.Menu(root)
    root["menu"] = menubar
    
    # Step 2 Make the pull-down menus on the menu bar
    
    
    change_menu = tkinter.Menu(menubar)
    menubar.add_cascade(menu = change_menu, label = "ChangeIt")
    
    show_menu = tkinter.Menu(menubar)
    menubar.add_cascade(menu= show_menu,label = "ShowIt")
    
    
    
    ## Step 3 Make menu items for each menu on the menu bar
    # Bind callbacks using lambda, as we have seen elsewhere,
    # but this time with a command= .. optional argument supplied.
    
    change_menu.add_command(label = "Increase the number",\
                            command = lambda:increase_number(data,1))


    change_menu.add_command(label = "Decrease the number",\
                            command = lambda:increase_number(data,-1))   
        
    show_menu.add_command(label = "SHow it in red",\
                          command = lambda: show(data,"red" ))
        
    show_menu.add_command(label = "Show it in yellow",\
                          command = lambda : show(data,"yellow"))
        
    root.mainloop()
    
def increase_number(data, amount):
    """
    Increase the number in the given Data object by the given amount and
    updates the Label in the given Data object that displays the Number

    """
    data.number = data.number +amount
    new_text  = "The number is {}".format(data.number)
    data.number_label["text"] = new_text
    
def show(data,color):
    new_text = "The number is {}".format(data.number)
    data.number_label["text"] = new_text
    data.number_label["background"] = color


main()


























    


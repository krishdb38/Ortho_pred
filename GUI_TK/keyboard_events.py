# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 08:45:57 2020

@author: krishna

"""


import tkinter
from tkinter import ttk

def main():
    root = tkinter.Tk()
    
    main_frame = ttk.Frame(root,padding = 20)
    main_frame.grid()
    
    left_button  = ttk.Button(main_frame,text = "left")
    left_button.grid()
    
    right_button = ttk.Button(main_frame, text = "right")
    right_button.grid()
    
    spin_button = ttk.Button(main_frame, text = "spin")
    spin_button.grid()
    
    left_button["command"] = lambda : go_left_button()
    right_button["command"]= lambda : go_right()
    spin_button["command"] = lambda : spin()
    
    
    root.bind_all("<KeyPress>", lambda event: pressed_a_key(event))
    root.bind_all("<KeyRelease>",lambda event : released_a_key(event))
    
    root.bind_all("<Key-L>", lambda event:go_left(event))
    root.bind_all("<Key-R>", lambda event: go_right(event))
    root.bind_all("<Key-r>", lambda event: go_right(event))
    root.bind_all("<Key-space>", lambda event: spin(event))

    root.mainloop()

def pressed_a_key(event):
    print("You pressed the", event.keysym, "key")
    

def released_a_key(event):
    print("you released the ",event.keysym, "key")
    
def go_left(event):
    print("You pressed the"+ event.keysym + "key:", end = "")
    print("go left")
    
def go_left_button():
    print("You clicked the Left butyton:")
    print("Go left")
    
def go_right(event = None):
    if event is None:
        print("Button press",end = "")
    else:
        print("You pressed teh "+ event.keysym+"key:",end = "")
    print("Go right")

def spin(event = None):
    if event is None:
        print("Button Press:",end = "")
    else:
        print("You pressed the "+event.keysym+ "key:",end = "")
    print("Spin")
    
main()
    
    















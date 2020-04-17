"""

Program developed by Krishna

"""

import tkinter
from tkinter import ttk
import time

def main():
    # Thus useal root and main Frame
    root = tkinter.Tk()
    mainframe = ttk.Frame(root, padding=20)
    mainframe.grid()
    
    # Checkbutton's and Radiobuttion's have their own labels
    checkbutton  = ttk.Checkbutton(mainframe,text = "Robuts rule")
    
    
    # Radiobutton's . We often put them onto a sub-frame
    # to group them visually. The "value" identifies which is selected
    
    radio_frame = ttk.Frame(mainframe, borderwidth = 10, relief = "groove")
    
    radio1 = ttk.Radiobutton(radio_frame, text = "Peter Pevensie,
                             value = "Peter")
    radio2 = ttk.Radiobutton(radio_frame, text = "Susan Pevensie", value = "sussan")
    
    radio3 = ttk.Radiobutton(radio_frame, text = "Edmund Pevensie",
                             value = "edmund")
    radio4 = ttk.Radiobutton(radio_frame, text = "Lucy Pevensie", value = "lucy")
    
    # This Button will show how it can interact with other widgets.
    
    # Checkbutton's and Radiobutton's can have an "observer" variable
    # that is bound to the state of the Checkbutton / Radiobutton.
    
    checkbutton_observer = tkinter.StringVar()
    checkbutton["variable"] = checkbutton_observer
    
    radio_observer  = tkinter.StringVar()
    
    for radio in  [radio1, radio2, radio3, radio4]:
        radio["variable"] = radio_observer # they all need the same observer
        
    # Bing callbacks using "command" and lambda, as we have seen elsewhere.
    checkbutton["command"] = lambda
    
    
    
    
    
    
main()
    
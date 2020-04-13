
import random

import tkinter as tk
import time
from time import strftime

bigLabels = ["Ortho_Detection"]
medLabels = ["Mode"]

window = tk.Tk()
#window.geometry("1600x800+0+0")
window.title("owPReMark_Orthologs_Detection_Software")

# Window logo in the top
photo = tk.PhotoImage(file = "./GUI_TK/icons/python2.png")
window.tk.call("wm","iconphoto",window._w,photo)
window.title("owPReMark, Software to detect Orthologes among Genomes")

text_input = tk.StringVar()
mode_input = tk.StringVar()

Tops = tk.Frame(window, width=1600,height =50,bg="powder blue")
Tops.pack(side = "top")

# frame 1 for main Software
f1 = tk.Frame(window,width = 800,height = 700,bg = "light blue").pack(side = "left")

# Frame 2 for extra Features
f2 = tk.Frame(window,width = 300,height = 700,bg= "light green").pack(side= "right")

prog_label = tk.Label(Tops,font = ("arial",45,"bold"),text = "owPReMark V1.0",fg = "blue")
prog_label.grid(row=0, column=0)

# !! Current Time
def time_():
    pass
    #string_ = strftime("%H:%M:%S %p")
    #time_label.config(text = string_)
    #time_label.after(1000, time)

local_time=time.asctime(time.localtime(time.time()))
time_label = tk.Label(Tops,font = ("arial",18,"bold"),text ="Started_Time-->"+local_time,fg = "white",bg = "purple", bd =10)
time_label.grid(row=1,column=0)

# time() for live time
window.mainloop()
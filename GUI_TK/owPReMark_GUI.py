
import random
import time
import tkinter as tk


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

# frame 1
f1 = tk.Frame(window,width = 800,height = 700,bg = "light blue").pack(side = "left")

# Frame 2
f2 = tk.Frame(window,width = 300,height = 700,bg= "light green").pack(side= "right")
# !! Current Time
local_time = time.asctime(time.localtime(time.time()))
prog_name = tk.Label(Tops,font = ("arial",40,"bold"),text = "owPReMark V1.0",fg = "blue", bd = 10)

prog_name.grid(row=0, column=0)



window.mainloop()
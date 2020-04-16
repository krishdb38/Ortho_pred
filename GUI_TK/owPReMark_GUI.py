
import random
import webbrowser
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import time
from time import strftime
import os


## !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
window = tk.Tk()
window.title("owPReMark_Orthologs_Detection_Software")
window.geometry("1000x800")
window.resizable(0,0)

photo = tk.PhotoImage(file = "./icons/python2.png")
#window.tk.call("wm","iconphoto",window._w,photo)
#window.title("owPReMark, Software to detect Orthologes among Genomes")

# !! ______________Defined_global_variables_Here_______________
## Lines
Tops = tk.Frame(window, width=1000,height =10,bg="green").pack(side="top")
buttoms = tk.Frame(window,width=1000,height =10,bg="green").pack(side="bottom")
# frame 1 for main Software
f1 = tk.Frame(window,width = 800,height = 1000,bg = "yellow")
# buttom_frame = tk.Frame(f1,width=800,height = 10,bg= "red").pack(side = "bottom")
f1.pack(side = "bottom")

# !!!  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Frame 2 for extra Features
f2 = tk.Frame(window,width = 500,height = 200,bg= "red",bd =15,).place(x=490,y=510)
f3 = tk.Frame(window,width = 500,height = 350,bg= "green").place(x=490,y=180)

# !!!!  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
prog_label = tk.Label(window,font = ("arial",45,"bold"),text = "owPReMark V1.0",fg = "blue")
prog_label.pack()

# !! Current Time
local_time=time.asctime(time.localtime(time.time()))
time_label = tk.Label(window,font = ("arial",18,"bold"),text ="Started"+local_time,fg = "white",bg = "purple", bd =10)
time_label.pack()

# !_______________________help & Credit Button__________________
def openURL(x):
    new = 2
    webbrowser.open(x,new = new)

def info():
    tk.messagebox.showinfo("info","Theragen Genome Care, South Korea,\n\
        Developer: Adhikari Krishna,\n\
        concept developed by : Kim\n\
        Supporting okg, bak")

tk.Button(window,text="Help....", width=10,padx=5,pady=8,bd=8,fg="green",font=('arial',16,"bold"),bg="powder blue",\
    command= lambda : openURL("https://www.google.com/")).place(x =800,y = 20)

tk.Button(window,text = "Credit", width=10,padx=1,pady=8,bd=8,fg="blue",font=('arial',16,"bold"),bg="powder blue",command = info).place(x =20,y = 20)
# Window logo in the top
# ###########################################################################
## !!__________________________Button_Program_run_________

def quit():
    msg = tk.messagebox.askquestion("Exit Application","Are you sure to exit the application",icon ="warning")
    if msg == "yes":
        window.destroy()
    else:
        tk.messagebox.showinfo("Return","you will now return to the application screen")

def reset():
    msg = tk.messagebox.askquestion("Are You sure to restore defaults","Are you sure to exit the application",icon ="warning")
    if msg == "yes":
        mode.set("0")
        inflation_factor.set("")
        select_folder.set("")
        cluster_result.set("")
    else:
        tk.messagebox.showinfo("Return","you will now return to the application screen")


def run_blast():
    pass

def mcl_clustering():
    pass

## !!XXXXXXXXXXX__Global_Variables__XXXXXXXXXXXXXXXXXXXXXXX
# Global variables
cpu_count = tk.IntVar()
blastp = tk.StringVar()
matrix = tk.StringVar()
blastp_data = tk.StringVar()
score_file = tk.StringVar()
threshold_score = tk.StringVar()

inflation_factor = tk.IntVar()
cluster_out = tk.StringVar()
USER_SELECTED_NUMBER= tk.IntVar() #later used split Function


# !_____________Create 3 button inside Fram1_

btn_runblast = tk.Button(f1,padx=16,pady=8,bd=16,fg="green",font=('arial',14,"bold"),width=10,text="Blastp_Only(1)",bg="powder blue")
btn_runblast.pack(side = "left")

btn_pre_blast = tk.Button(f1,padx=16,pady=8,bd=16,fg="green",font=('arial',14,"bold"),width=12,text="Check_blastp(2)",bg="powder blue")
btn_pre_blast.pack(side = "left")

btn_pre_blast = tk.Button(f1,padx=16,pady=8,bd=16,fg="green",font=('arial',14,"bold"),width=12,text="MCL_Clustering(3)",bg="powder blue")
btn_pre_blast.pack(side = "left")

btnreset=tk.Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',14,"bold"),width=8,text="Reset **",bg="powder blue",command=reset)
btnreset.pack(side = "left")

btnexit=tk.Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',14,"bold"),width=8,text="Exit !!",bg="powder blue",command=quit)
btnexit.pack(side = "left")
## !!____________________________________labels & Entry_______________________________________________

inflation_label = tk.Label(window,text = "Threshold Value",font=('arial',16,"bold")).place(x=20,y=200)
inflation_entry = tk.Entry(window,bd =5,).place(x=300,y=200)

inflation_label = tk.Label(window,text = "Enter inflation factor",font=('arial',16,"bold")).place(x=20,y=250)
inflation_entry = tk.Entry(window,bd =5,).place(x=300,y=250)

clustering_label = tk.Label(window,text = "Name for Clustering File",font=('arial',16,"bold")).place(x=20,y=300)
clustering_entry = tk.Entry(window,bd =5,).place(x=300,y=300)

clustering_label = tk.Label(window,text = "Name for Log File",font=('arial',16,"bold")).place(x=20,y=350)
clustering_entry = tk.Entry(window,bd =5).place(x=300,y=350)

cpu_label = tk.Label(window, text = "Number of Cpu.",font=('arial',16,"bold")).place(x=20,y=400)
tk.Label(window,text = "less than"+str(os.cpu_count()),font=('arial',12,"bold"),fg="red").place(x=200,y=400)

clustering_entry = tk.Entry(window,bd =5).place(x=300,y=400)


# !!__________________________Browse_Button__________________________________________
species_folder = tk.StringVar()
blastp_folder = tk.StringVar()

def browse_button(arg):
    "Allow user to select a directory and store it in global variable variable_name --folder_path"
    global species_folder, blastp_folder 
    filename = filedialog.askdirectory()
    arg.set(filename)
    return arg

tk.Label(window,textvariable = species_folder).place(x=100,y=500)
tk.Button(text="Species",command =lambda: browse_button(species_folder)).place(x=20,y=500)


tk.Label(window,textvariable = blastp_folder).place(x=100,y=550)
tk.Button(text="blastp_file",command =lambda: browse_button(blastp_folder)).place(x=20,y=550)



## !! _____________________Display input_parameter inside frame 3 with button
def display(input):
    pass

show_variable = tk.Button(f3,text = "show_variable").pack(side = "left")


# !!___________________________________Radio_Button_for Matrix_selection_________________________________
v = tk.StringVar()

def show_choice():
    print(v.get())
tk.Label(window,text = "Choose a Metrix for blastp",font=('arial',16,"bold")).place(x=50,y=600)

tk.Radiobutton(window,text = "BLOSUM45",variable=v,value = "BLOSUM45",command = show_choice,font=('arial',10,"bold")).place(x=10,y=650)
tk.Radiobutton(window,text = "BLOSUM62",variable = v,value="BLOSUM62",command = show_choice,font=('arial',10,"bold")).place(x=120,y=650)
tk.Radiobutton(window,text = "BLOSUM82",variable = v,value="BLOSUM82",command = show_choice,font=('arial',10,"bold")).place(x=240,y=650)

window.mainloop()

if __name__ == '__main__':
    pass
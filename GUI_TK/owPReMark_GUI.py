
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
window.geometry("1000x800+100+20")
window.resizable(0,0)

photo = tk.PhotoImage(file = "./icons/python2.png")
window.tk.call("wm","iconphoto",window._w,photo)
window.title("owPReMark, Software to detect Orthologes among Genomes")

# !! ______________Defined_global_variables_Here_______________
## Lines
Tops = tk.Frame(window, width=1000,height =10,bg="green").pack(side="top")
buttoms = tk.Frame(window,width=1000,height =10,bg="green").pack(side="bottom")
# frame 1 for main Software
f1 = tk.Frame(window,width = 800,height = 100,bg = "yellow")
# buttom_frame = tk.Frame(f1,width=800,height = 10,bg= "red").pack(side = "bottom")
f1.pack(side = "bottom")

# !!!  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Frame 2 for extra Features
f2 = tk.Frame(window,width = 500,height = 200,bg= "#FEF8DD",bd =15,).place(x=490,y=510)
f3 = tk.Frame(window,width = 500,height = 350,bg= "#8DE4FD").place(x=490,y=180)

# !!!!  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
prog_label = tk.Label(window,font = ("arial",40,"bold"),text = "owPReMark V1.0",fg = "white",bg="#093145")
prog_label.pack()

# !! Current Time
local_time=time.asctime(time.localtime(time.time()))
time_label = tk.Label(window,font = ("arial",18,"bold"),\
                      text ="Started-->  "+local_time,fg = "#093145", bd =10)
time_label.pack()

local_time=time.asctime(time.localtime(time.time()))
time_label = tk.Label(window,font = ("arial",18,"bold"),\
                      text ="Consumed-->  "+local_time,fg = "#093145", bd =10)
time_label.pack()

# !_______________________help & Credit Button__________________
def openURL(x):
    """Open the help menu of Git Hub."""
    new = 2
    webbrowser.open(x,new = new)

def info():
    """Show the information about company & developer."""
    tk.messagebox.showinfo("info","Theragen Genome Care, South Korea,\n\
                           concept developed by : Kim\n\
                        Developer: Adhikari Krishna,\n\
                        Supporting bak, okg")

tk.Button(window,text="Help....", width=10,padx=5,pady=8,bd=8,font=('arial',16,"bold"),bg="powder blue",\
    command= lambda : openURL("https://www.google.com/")).place(x =800,y = 20)

tk.Button(window,text = "Credit", width=10,padx=1,pady=8,bd=8,font=('arial',16,"bold"),bg="powder blue",command = info).place(x =20,y = 20)
# Window logo in the top
# ###########################################################################
## !!__________________________Button_Program_run_________

def quit():
    """Exit the program."""
    msg = tk.messagebox.askquestion("Exit Application","Are you sure to exit the application",icon ="warning")
    if msg == "yes":
        window.destroy()
    else:
        tk.messagebox.showinfo("Return","you will now return to the application screen")

def reset():
    """Reset the defaults parameters."""
    msg = tk.messagebox.askquestion("Are You sure to restore defaults","Are you sure to exit the application",icon ="warning")
    if msg == "yes":

        blastp_data. set("./blastp_data/")
        score_file.set("./score_file/")
        threshold_score.set(5)
        inflation_factor.set(1.5)
        log_file.set("log")
        cluster_out.set("results")
        cpu_count.set(1)
        
    else:
        tk.messagebox.showinfo("Return","you will now return to the application screen")


def run_blast():
    """Run blastp in the system."""
    pass

def mcl_clustering():
    """Perform MCL_clustering."""
    pass

## !!XXXXXXXXXXX__Global_Variables__XXXXXXXXXXXXXXXXXXXXXXX
# Global variables
cpu_count = tk.IntVar()
cpu_count.set(1)

blastp_data = tk.StringVar()
blastp_data. set("./blastp_data/")

score_file = tk.StringVar()
score_file.set("./score_file/")

threshold_score = tk.IntVar()
threshold_score.set(5)

inflation_factor = tk.IntVar()
inflation_factor.set(1.5)

log_file = tk.StringVar()
log_file.set("log")

cluster_out = tk.StringVar()
cluster_out.set("results")

matrix_selection = tk.StringVar()

# Setting the default value



# !_____________Create 3 button inside Fram1_

btn_runblast = tk.Button(f1,padx=16,pady=8,bd=16,fg="green",\
                         font=('arial',14,"bold"),width=10,text="Blastp_Only(1)",\
                             bg="powder blue").pack(side = "left")


btn_pre_blast = tk.Button(f1,padx=16,pady=8,bd=16,fg="green",\
                          font=('arial',14,"bold"),width=12,text="Check_blastp(2)",\
                              bg="powder blue").pack(side = "left")


btn_pre_blast = tk.Button(f1,padx=16,pady=8,bd=16,fg="green",\
                          font=('arial',14,"bold"),width=12,\
                              text="MCL_Clustering(3)",\
                                  bg="powder blue").pack(side = "left")


btnreset=tk.Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',14,"bold"),\
                   width=8,text="Reset **",bg="powder blue",\
                       command=reset).pack(side = "left")

    
btnexit=tk.Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',14,"bold"),width=8,text="Exit !!",bg="powder blue",command=quit)
btnexit.pack(side = "left")
## !!____________________________________labels & Entry_______________________________________________

threshold_label = tk.Label(window,text = "Threshold Score",font=('arial',16,"bold")).place(x=20,y=200)
threshold_entry = tk.Entry(window,bd =5,bg = "#CAF1DE", textvariable = threshold_score,font=('arial',12,"bold"), width=16).place(x=300,y=200,height = 35)

inflation_label = tk.Label(window,text = "Enter inflation factor",font=('arial',16,"bold")).place(x=20,y=250)
inflation_entry = tk.Entry(window,bd =5,bg = "#CAF1DE",textvariable = inflation_factor,font=('arial',12,"bold"), width=16).place(x=300,y=250,height = 35)

clustering_label = tk.Label(window,text = "Name for Clustering File",font=('arial',16,"bold")).place(x=20,y=300)
clustering_entry = tk.Entry(window,bd =5,bg = "#CAF1DE",textvariable = cluster_out,font=('arial',12,"bold"), width=16).place(x=300,y=300,height = 35)

log_label = tk.Label(window,text = "Name for Log File",font=('arial',16,"bold")).place(x=20,y=350)
log_entry = tk.Entry(window,bd =5,bg = "#CAF1DE", textvariable = log_file,font=('arial',12,"bold"), width=16).place(x=300,y=350,height = 35)

cpu_label = tk.Label(window, text = "Number of Cpu.",font=('arial',16,"bold")).place(x=20,y=400)
tk.Label(window,text = "less than"+str(os.cpu_count()),font=('arial',12,"bold"),fg="red").place(x=200,y=400)

cpu_entry = tk.Entry(window,bd =5, textvariable =cpu_count,bg = "#CAF1DE", font=('arial',12,"bold"), width=16).place(x=300,y=400,height = 35)


# !!__________________________Browse_Button__________________________________________
species_folder = tk.StringVar()
blastp_folder = tk.StringVar()

def browse_button(arg):
    """Allow user to select a directory and store it in global variable variable_name --folder_path."""
    global species_folder, blastp_folder 
    filename = filedialog.askdirectory()
    arg.set(filename)
    return arg

tk.Label(window,textvariable = species_folder).place(x=100,y=500)
tk.Button(window,text="Species_folder",font = ("arial",12,"bold"),bd =2,
          bg ="#8DE4FD",command =lambda: browse_button(species_folder)).place(x=20,y=500)


tk.Label(window,textvariable = blastp_folder).place(x=100,y=550)

tk.Button(window,text="blastp_folder",font = ("arial",12,"bold"),bd =2,
          bg ="#8DE4FD",command =lambda: browse_button(blastp_folder)).place(x=20,y=550)


## !! _____________________Display input_parameter inside frame 3 with button
def display():
    """Display the information about variable."""
    
    print(score_file.get())
    info_lbl.config(text = "this is test")
    print(inflation_factor.get())
    print(log_file.get())
    print(cpu_count.get())
    
info_lbl = tk.Label(f2).place(x=490,y=540)

show_variable = tk.Button(f2,text = "show_IN_variable",command = display).place(x=490,y=510)
#info_lbl = tk.Label(f2,text="Selected Variables are:").place(x=490,y=540)



# !!___________________________________Radio_Button_for Matrix_selection_________________________________
v = tk.StringVar()

def show_choice():
    
    global matrix_selection
    matrix_selection = v.get()
    print(v.get())
    
tk.Label(window,text = "Choose a Matrix for blastp",font=('arial',14,"bold")).place(x=50,y=620)

tk.Radiobutton(window,text = "BLOSUM45",variable=v,value = "BLOSUM45", \
               command = show_choice,font=('arial',10,"bold")).place(x=10,y=650)
tk.Radiobutton(window,text = "BLOSUM62",variable = v,value="BLOSUM62",\
               command = show_choice,font=('arial',10,"bold")).place(x=120,y=650)
tk.Radiobutton(window,text = "BLOSUM82",variable = v,value="BLOSUM82",command = show_choice,font=('arial',10,"bold")).place(x=240,y=650)

window.mainloop()

if __name__ == '__main__':
    pass
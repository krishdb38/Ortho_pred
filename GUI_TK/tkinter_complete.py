import tkinter as tk
from tkinter import messagebox


window = tk.Tk()

# lets fix the windows  size
#window.geometry("500x500")

# Window logo in the top
#photo = tk.PhotoImage(file = "./GUI_TK/icons/python2.png")
#window.tk.call("wm","iconphoto",window._w,photo)
window.title("owPReMark, Software to detect Orthologes among Genomes")

# !! ____________________Buttot_______________________________
def hello_callback():
    return messagebox.showinfo("Hello Python","Wrong_ Input")

B = tk.Button(window, text= "hello", command = hello_callback,)
B.place(x = 50,y =50)
# !!_______________________Canvas__________________________________
# Inside canvas image , line , oval , polygon
# C = tk.Canvas(window,bg = "blue",height = 250)
# coord = 10,50,240,210
# arc = C.create_arc(coord,start = 0, extent = 150 , fill = "red")
# line = C.create_line(10,10,200,200,fill = "white")
# C.pack()
 # !!_______________________Check_Buttotn__________________________
# window = tk.Checkbutton(master, options,....)
# options are activebackground, activeforeground, bg, bitmap, command,
# cursor, disabledforeground, font, fg, height, highlightcolor, image,
# justify ,  offvalue, offvalue, onvalue, padx, pady, relief, selectcolor,
# selectimage, state, text, underline, variable , width, wraplength,
# Methods-- deselect(), flash(), invoke(), select(), toggle()
CheckVar1 = tk.IntVar()
CheckVar2 = tk.IntVar()
C1 = tk.Checkbutton(window, text ="Blastp_New",variable = CheckVar1,\
    onvalue =1, offvalue = 0, height= 5,\
        width =20,)

C2 = tk.Checkbutton(window, text = "Blastp_precalculated",variable =CheckVar2,\
    onvalue= 1, offvalue = 0, height = 5,width=20)
C1.pack()
C2.pack()
# !!____________________________Entry________________________________
# Entry widget is used to display a single-line text field
# for accepting values from user
# window = tk.Entry(master,option)
# options are -- bg, bd, command, cursor, font, exportselection, fg, highlightcolor,
# justify relief selectbackground selectborderwidth selectforeground show ,
#  state, textvariable, width, xscrollcommand
# Methods delete(first , last =None), get(), icursor(index),insert(index, s),
# select_adjust(index), select_clear(), select_from(index), select_present()
# select_range(start, end), select_to(index), xview(index), xview_scroll(number, what)
L1 = tk.Label(window, text = "Developed_by_Krish")
L1.pack(side = "left")
E1 = tk.Entry(window, bd = 5)
E1.pack(side = "left")

# !!_______________________Frame___________________________________________
# Frame widget is very important for the process of grouping and organizing
# other widgets in a some how friendly wat.works like a container
# w = tk.Frame(master, options.....)
# options are bg, bd, cursor, height, highlightbackground, highlightcolor, 
# highlightcolor, highlightthickness, relief, width

frame = tk.Frame(window).pack()
bottomframe = tk.Frame(window).pack(side = "left")

redbutton = tk.Button(frame, text = "Red", fg = "black").pack(side = "left")
greenbutton = tk.Button(frame, text = "Brown", fg = "brown").pack(side = "left")

bluebutton = tk.Button(frame, text = "Blue", fg = "blue").pack(side = "left")

blackbutton = tk.Button(bottomframe,text = "Black", fg = "black").pack(side = "bottom")
# !______________________Label_______________________________________________
# w = Label(master, option)
# options--> anchor, bg, bitmatp, bd, cursor, font, fg , height, image, 
# justify, padx, pady , relief, text, textvariable, underline, width,
# wraplength, 

var3 = tk.StringVar()
label = tk.Label(frame,textvariable = var3, relief = "raised")
var3.set("This program will Find orthologs")
label.pack(side = "left")

# !_________________________Listbox_____________________
# is used to display a list of items from which a user can select a number of items
# w = tk.Listbox(master, options...)
# options --> bg, bd, cursor, font, fg, height,\
# highlightcolor, highlightthickness, relief, selectbackground, selectmode
# width, xscrollcommand, yscrollcommand,
# Methods are --> activate(index),curselection(), delete(first, last = None)
# get(first, last= None), index(i), insert(index, *elements), nearest(y),
# see(index), size(), xview(), xview_moveto(fraction)
# xview_scroll( number, what), yview(), yview_moveto(fraction)
# yview_scroll(number, what)
Lb1 = tk.Listbox(frame)
Lb1.insert(1,"Python")
Lb1.insert(2,"Perl")
Lb1.insert(3,"C")
Lb1.insert(4,"PHP")
Lb1.insert(5,"JSP")
Lb1.insert(6,"Rupy")
Lb1.pack()

# !___________________Menubutton_______________________
# tk.Menubutton(master, options,......)
# options:--> activebackground, activeforeground, anchor, bg, bitmap, bd, cursor, direction,
# disabledforeground, fg, height, highlightcolor, image, justify, menu, padx, pady, relief
# state, text, textvariable, underline, width, wraplength

mb = tk.Menubutton(window, text= "condiments", relief = "raised")
# !! mb.grid()
mb.menu = tk.Menu(mb,tearoff = 0)
mb["menu"] = mb.menu

mayoVar = tk.IntVar()
ketchVar = tk.IntVar()

mb.menu.add_checkbutton(label = "mayo",\
    variable= mayoVar)
mb.menu.add_checkbutton( label = "ketchup",\
    variable = ketchVar)
mb.pack()

## !!!______________________Menu_____________________
# The goal of this widget is to allow us to create all kinds of menus
# that can be used by our app.3 types pop-up, toplevel and pull-down
# w = tk.Menu(master, options)
# options are :--> activebackground, activeborderwidth, activeforeground,
# bg, bd, cursor, disabledforeground, font, fg, postcommand, relief,
# image, selectcolor, tearoff, title
# methods add_command(options), add_radiobutton(options),add_checkbutton(options)
# add_cascade(options), add_separator(), add(type, options), delete(startindex[,endindex])
# entryconfig(index, options), index(item), insert_separator(index),invoke(index)
# type (index)
def donothing():
    filewin = tk.Toplevel(window)
    button = tk.Button(filewin,text= "Do nothing button")
    button.pack()


menubar = tk.Menu(window)

filemenu = tk.Menu(menubar,tearoff = 0)
filemenu.add_command(label = "New",command = donothing)
filemenu.add_command(label = "open",command = donothing)
filemenu.add_command(label = "save",command = donothing)
filemenu.add_command(label = "Save as..",command = donothing)
filemenu.add_separator()
filemenu.add_command(label = "Exit",command = window.quit)
menubar.add_cascade(label = "File",menu = filemenu)

editmenu = tk.Menu(menubar, tearoff= 0)
editmenu.add_command(label = "Undo", command = donothing)
editmenu.add_separator()


editmenu.add_command(label = "Cut", command = donothing)
editmenu.add_command(label = "Copy", command = donothing)
editmenu.add_command(label = "Paste", command = donothing)
editmenu.add_command(label = "Delete", command = donothing)
editmenu.add_command(label = "Select All", command = donothing)

#menubar.add_cascasde(label = "Edit", menu = editmenu)
helpmenu = tk.Menu(menubar, tearoff =0)
helpmenu.add_command(label= "help Index", command = donothing)
helpmenu.add_command(label= "Help",)
 # !____________________Message_______________________________
# w = tk.Message(master, options)
# option & Descriptions
# anchor, bg, biotmap, bd, cursor, font, fg, fheight, image, 
# justify, padx, pady, relief, text, textvariable,
# underline, width, wraplength
var = tk.StringVar()
label = tk.Message(window, textvariable = var, relief = "raised")
var.set("Hey !? How are doing")

label.pack()
# !________________________Radiobutton____________________________
# w = tk.Radiobutton(master, option,......)
# activebackground, activeforeground, anchor , bg, bitmap, borderwidth,command
# cursor, font, fg, height, highlightbackground, highlightcolor, image, justify
# padx, pady, relief, selectcolor, selectimage, state, text, textvariable,
# underlilne, value, variable, width, wraplength,
# Methods are --> deselect(),flash(), invoke(), select()
def sel():
    selection = "You selected the option",str(var.get())
    label.config(text = selection)

var = tk.IntVar()
R1 = tk.Radiobutton(window, text = "Blastp", variable = var,\
    value = 1,command = sel).pack(anchor = "w")
R2 = tk.Radiobutton(window, text = "Blastp with Precalculated",\
    variable = var, value =2,command = sel).pack(anchor = "w")
R3 = tk.Radiobutton(window, text = "Clustering with MCL",\
    variable = var, value = 3, command = sel).pack(anchor = "w")
label = tk.Label(window).pack()
# !_______________________Scale__________________________
# Provides a graphical slider object that allows you to select values from a specific scale
# w = tk.Scale(master, options)
#options: --> activebackground, bg, bd, command, cursor, digits, font, fg,
# from_ , highlightbackground, highlightcolor, label , length , orient , relief ,
# repeatdelay, resolution , showvalue, sliderlength , state, takefocus, tickinterval
# to , throughcolor, variable, width
# Methods
def sel():
    selection = "Value = "+ str(var.get())
    label.config(text = selection)

var = tk.DoubleVar()
scale = tk.Scale(window, variable = var,bg="red")
scale.pack(anchor = "center")
button = tk.Button(window, text = "Get Scale Value", command = sel)
button.pack(anchor = "center")

label = tk.Label(window).pack()
window.config(menu = menubar)

# !! _____________________Scrollbar____________________________
# a slide controller that is used to implement vertical scrolled widgets,
# such as Listbox, Text, and Canvas.Note, can also creatre horizontal scrtollbars on Entry widgets.
# w = Scrollbar(master, option,.....)
# options --> activebackfround, bg, bd, command, cursor, elementborderwidth,
# highlightbackground, highlightcolor, highlightthickness, jump
# orient, repeatdelay , repeatinterval, takefocus, troughcolor, width
# Methods --> get() , set(first, last)

# scrollbar = tk.Scrollbar(window)
# scrollbar.pack(side = "right", fill = "y")

# mylist = tk.Listbox(window, yscrollcommand = scrollbar.set)
# for line in range(5):
#     mylist.insert("end","This is line number"+str(line))

# mylist.pack(side = "left", fill = "both")
# scrollbar.config(command = mylist.yview)

# !______________________text__________________________________
# w = tk.Text(master, options.....)
# options and descriptions
# bg, bd, cursor, exportselection, font, fg, height, highlightbackground,
# highlightcolor, hightlightthickness, insertbackground, insertborderwidth,
# insertofftime, insertontime, insertwidth, padxx, pady, relief, selectbackground,
# selectborderwidth, spacing1, spacing2, spacing3, state, tabs, width, wrap,
# xscrollcommand, yscrollcommand,
# Methods & Description 
#delete(startindex[,endindex]), get(startindex[,endindex]), index(index)
# insert(indwex[,string]..), see(index)

# https://www.tutorialspoint.com/python3/tk_text.htm

# text = tk.Text(window)
# text.insert("insert", "Hello....")
# text.insert("end","BYE Bye Bye....")
# text.pack()

# !___________________Toplevel____________________
# toplevel widgets work as windows that are directly managed by the window manager
# They do not necessarily have  a parent widget on top of them
# w = tk.Toplevel(options, ....)
# parameters --> bg, bd, cursor, class_ , font, fg, height, relief, width
# Methods --> deiconify(), frame() , group(window), iconify(), protocol(name, function)
# iconify(), state(), transient([master]), withdraw(), maxsize(width, height)
# state(), transient([master]), withdraw(), maxsize(width, height),minsize(width, height)
# positionfrom(who), resizable(width, height), sizefrom(who), title(string)

# top = tk.Toplevel()
# top.title("This is title of top")
# top.mainloop()

# !!_________________Spinbox_______________________________
# w = tk.Spinbox(master, option,.....)
# option & Descriptions
# activebackground, bg, bd, command, cursor, disablebackground, disabledforeground, fg
# font, format, from_, justify, relief, repeatdelay, repeatinterval, state,
# textvariable, to, validate, validatecommand, values, vcmd, wrap, xscrollcommand
# Methods --> delete(startindex[,endindex]) , get(startindex[,endindex]), identify
# index(index), insert(index[,string]), invoke(element), 
w = tk.Spinbox(window,from_ = 0, to = 10).pack(side = "top")


# !!___________________Paned_WIndow_________________________________
# PanedWindow is a container widget that may contain any number of paes, arranged horizontally or vertically.
# w = tk.PanedWindow(master, options....)
# Parameters bg, bd, borderwidth, cursor, handlepad, height, orient, relief, sashcurser, sashrelief, sashwidth,
# showhandle, width. Methods --> add(child, options), get(startindex[,endindex]), config(options)
m1 = tk.PanedWindow()
m1.pack(fill = "both", expand = 1)

left = tk.Entry(m1,bd = 5)
m1.add(left)

m2 = tk.PanedWindow(m1, orient = "vertical")
m1.add(m2)

top = tk.Scale(m2, orient = "horizontal")
m2.add(left)
# !!__________________LabelFGrame______________
# labelframe is a simple container widget.
# Its primary purpose is to act as a spacer or container for complex window layouts
# w = tk.LabelFrame(master, option,.......)
# options are bg, bd, cursor, font, height, labelAnchor, highlightbackground, highlightcolor,
# highlightthickness, relief, text, width, 
label_frame = tk.LabelFrame(window, text = "This is a labelFrame")
label_frame.pack(fill= "both",expand = "yes")

left = tk.Label(label_frame, text = "Inside the label Frame")
left.pack()

# !!____________________tkMessageBox__________
# is used to display message boxes in your applications.
# This module provides a number of functions that you can use to display an appropriate message.
# tkMessageBox.FunctionName(title, message[,options])
# parameters ***
# FunctionName --is the name of the appropriate message box function
# title --> text to be displayed in the title bar of a message box.
# message --> is the text to be displayed as a message
# options --> are alternative choices that you may use to tailor a standard message box.
# can use one of the following functions with dialogue box
# showinfo(), showwarning(), showerror(), askquestion(), askokcancel(), askyesno() , askretrycancel()
def hello():
    tk.messagebox.showinfo("Say__Hello","Hello World")
    tk.messagebox.showwarning("Not_match","Check_the_File")
B1 = tk.Button(window, text = "Say Hello", command = hello)
B1.place(x =0,y=0)




window = tk.Checkbutton(window,)
window.mainloop()



































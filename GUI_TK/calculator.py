import tkinter as tk

root = tk.Tk()
values = "0"

temp = ""
inputs = tk.StringVar()
inputs.set(values)

def dot(num):
    global values
    count = 0
    for i in range(len(values)):
        if values[-i-1] in "+-*/":
            count = i
            break
        if "." in values[-count]:
            pass
        else:
            values = values + num
        inputs.set(values)
        
        
        
def zero(any):
    global values
    if(values)=="0":
        pass
    if values[-1] in "+-*/" :
        values = values + any
        inputs.set(values)
    elif len(values) == 1 and values[0] =="1":
        values = values +any
        inputs.set(values)
        pass
    elif len(values)>1:
        if (values[-2] in "+-/*" and values[-1] =="0"):
            pass
        else:
            values = values + any
            inputs.set(values)
            
    else:
        pass
    
def inputValue(any):
    global values
    if values[0] =="0":
        values = str(any)
    elif len(values) >1 and (values[-2] in "+-*/" and values[-1] =="0"):
        values = values[:-1]+str(any)
        
    else:
        values = values +str(any)
    
    inputs.set(values)
    
def inputSign(any):
    global values 
    if values =="0" and any == "-":
        values = str(any)
        inputs.set(values)
        
    if values[-1] in "./+-*":
        values = values[:-1]+str(any)
        inputs.set(values)
        
    else:
        values = values +str(any)
        inputs.set(values)
        
        
def delete(any):
    global values
    values = values[:-1] # delete by 1 position according to last index
    inputs.set(values)
    
# Restore the input
def clear(any):
    global values
    values = "0"
    inputs.set(values)
    
def calculate(any):
    try:
        global values
        result = str(eval(values))
        inputs.set(result)
        values = result
        
    except:
        inputs.set("Error")
        values = ""

        
        
### Develope the GUI 
root.geometry("300x320")
root.resizable(0,0) #Fixed the window
display = tk.Entry(root, textvariable = inputs, font = ("arial", 15,"bold"), justify = "right", width =15)

display.grid(row =0,column =1, columnspan =4, ipadx = 15, ipady = 10)

tk.Label(root, text ="krishna").grid(row = 5,column=0)
# ! **************

buttonBack = tk.Button(root,text = "Back",command = lambda: delete("<-"),width = 5, height = 2)
buttonBack.grid(row =1 ,column=1) 

 
root.mainloop()




















































        
        
        
        
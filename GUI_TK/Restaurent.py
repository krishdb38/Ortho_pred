# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 14:39:38 2018

@author: Deepak
"""

from tkinter import *
import random
import time;

window=Tk()
window.geometry("1600x800+0+0")
window.title("Restaurent Management System")

text_input=StringVar()
operator=""

Tops=Frame(window,width=1600,height=50,bg="powder blue")
Tops.pack(side=TOP)


f1=Frame(window,width=800,height=700)
f1.pack(side=LEFT)


f2=Frame(window,width=300,height=700)
f2.pack(side=RIGHT)
'''==============================================Time============================================================'''
localtime=time.asctime(time.localtime(time.time()))
'''=============================================Time Info======================================================='''
lblInfo=Label(Tops,font=("arial",20,"bold"),text="owPReMark",fg="red",bd=10)
lblInfo.grid(row=0,column=0)

lblInfo=Label(Tops,font=("arial",20,"bold"),text=localtime,fg="red",bd=10)
lblInfo.grid(row=1,column=0)
'''============================================Calculator======================================================='''
def btnclick(numbers):
    global operator
    operator = operator + str(numbers)
    text_input.set(operator)
def btnclear():
    global operator
    operator=""
    text_input.set("")
def btnequal():
    global operator
    sumpup=str(eval(operator))
    text_input.set(sumpup)
    operator=""
def ref():
    x=random.randint(500,50000000)
    randstr=str(x)
    rand.set(randstr)
    
    cof=float(fries.get())
    cob=float(burger.get())
    cot=float(tikki.get())
    coc=float(chicken.get())
    coch=float(cheese.get())
    cod=float(drinks.get())
    cos=float(shakes.get())
    
    costoffries=cof*30
    costofburger=cob*40
    costoftikki=cot*20
    costofchicken=coc*80
    costofcheese=coch*50
    costofdrinks=cod*30
    costofshakes=cos*50
    
    CostofMeal="Rs",str("%.2f" % (costoffries+costofburger+costoftikki+costofchicken+costofcheese+costofdrinks+costofshakes))
    
    totalcost=costoffries+costofburger+costoftikki+costofchicken+costofcheese+costofdrinks+costofshakes
    Tax=((costoffries+costofburger+costoftikki+costofchicken+costofcheese+costofdrinks+costofshakes)*0.18)
    
    Costoftax="Rs",str("%.2f" % Tax)
    
    Total="Rs",str("%.2f" % (Tax+totalcost))
    meal.set(CostofMeal)
    gst.set(Costoftax)
    subtotal.set(CostofMeal)
    total.set(Total)
def qexit():
    window.destroy()
def reset():
    rand.set("")
    fries.set("")
    burger.set("")
    tikki.set("")
    chicken.set("")
    cheese.set("")
    drinks.set("")
    meal.set("")
    shakes.set("")
    gst.set("")
    subtotal.set("")
    total.set("")
txtDisplay=Entry(f2,font=('arial',20,'bold'),textvariable=text_input,bd=30,insertwidth=4,bg="powder blue",justify='right')
txtDisplay.grid(columnspan=4)
'''====================================================================================='''
btn7=Button(f2,text=7,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(7))
btn7.grid(row=2,column=0)

btn8=Button(f2,text=8,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(8))
btn8.grid(row=2,column=1)

btn9=Button(f2,text=9,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(9))
btn9.grid(row=2,column=2)

Addition=Button(f2,text='+',padx=14,pady=14,bd=8,fg="red2",font=('arial',20,"bold"),command=lambda: btnclick("+"))
Addition.grid(row=2,column=3)
'''========================================================================================='''
btn4=Button(f2,text=4,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(4))
btn4.grid(row=3,column=0)

btn5=Button(f2,text=5,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(5))
btn5.grid(row=3,column=1)

btn6=Button(f2,text=6,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(6))
btn6.grid(row=3,column=2)

Subtract=Button(f2,text='-',padx=14,pady=14,bd=8,fg="red2",font=('arial',20,"bold"),command=lambda: btnclick("-"))
Subtract.grid(row=3,column=3)
'''==========================================================================================='''
btn1=Button(f2,text=1,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(1))
btn1.grid(row=4,column=0)

btn2=Button(f2,text=2,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(2))
btn2.grid(row=4,column=1)

btn3=Button(f2,text=3,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(3))
btn3.grid(row=4,column=2)

Multiply=Button(f2,text='*',padx=14,pady=14,bd=8,fg="red2",font=('arial',20,"bold"),command=lambda: btnclick("*"))
Multiply.grid(row=4,column=3)
'''============================================================================================='''
btn0=Button(f2,text=0,padx=14,pady=14,bd=8,fg="black",font=('arial',20,"bold"),command=lambda: btnclick(0))
btn0.grid(row=5,column=0)

btnc=Button(f2,text="C",padx=14,pady=14,bd=8,fg="red2",font=('arial',20,"bold"),command=btnclear)
btnc.grid(row=5,column=1)

btn=Button(f2,text="=",padx=14,pady=14,bd=8,fg="red2",font=('arial',20,"bold"),command=btnequal)
btn.grid(row=5,column=2)
Divide=Button(f2,text='/',padx=14,pady=14,bd=8,fg="red2",font=('arial',20,"bold"),command=lambda: btnclick("/"))
Divide.grid(row=5,column=3)
'''=============================================Restaurent info 1============================================'''
rand=StringVar()
fries=StringVar()
burger=StringVar()
tikki=StringVar()
chicken=StringVar()
cheese=StringVar()
drinks=StringVar()
meal=StringVar()
shakes=StringVar()
gst=StringVar()
subtotal=StringVar()
total=StringVar()
lblrefrence=Label(f1,font=('arial',16,'bold'),text='Reference',bd=16,anchor='w')
lblrefrence.grid(row=0,column=0)

txtrefrence=Entry(f1,font=('arial',16,'bold'),textvariable=rand,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtrefrence.grid(row=0,column=1)

lblfries=Label(f1,font=('arial',16,'bold'),text='Fries',bd=16,anchor='w')
lblfries.grid(row=1,column=0)

txtfries=Entry(f1,font=('arial',16,'bold'),textvariable=fries,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtfries.grid(row=1,column=1)

lblburger=Label(f1,font=('arial',16,'bold'),text='Burger Meal',bd=16,anchor='w')
lblburger.grid(row=2,column=0)

txtburger=Entry(f1,font=('arial',16,'bold'),textvariable=burger,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtburger.grid(row=2,column=1)

lbltikki=Label(f1,font=('arial',16,'bold'),text='Allo Tikki ',bd=16,anchor='w')
lbltikki.grid(row=3,column=0)

txttikki=Entry(f1,font=('arial',16,'bold'),textvariable=tikki,bd=10,bg="powder blue",justify="right",insertwidth=4)
txttikki.grid(row=3,column=1)

lblchicken=Label(f1,font=('arial',16,'bold'),text='Chicken Meal',bd=16,anchor='w')
lblchicken.grid(row=4,column=0)

txtchicken=Entry(f1,font=('arial',16,'bold'),textvariable=chicken,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtchicken.grid(row=4,column=1)

lblcheese=Label(f1,font=('arial',16,'bold'),text='Cheese Meal',bd=16,anchor='w')
lblcheese.grid(row=5,column=0)

txtcheese=Entry(f1,font=('arial',16,'bold'),textvariable=cheese,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtcheese.grid(row=5,column=1)
'''==========================================Restaurent info 2==============================================='''
lbldrinks=Label(f1,font=('arial',16,'bold'),text='Drinks',bd=16,anchor='w')
lbldrinks.grid(row=0,column=2)

txtdrinks=Entry(f1,font=('arial',16,'bold'),textvariable=drinks,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtdrinks.grid(row=0,column=3)

lblShakes=Label(f1,font=('arial',16,'bold'),text='Shakes',bd=16,anchor='w')
lblShakes.grid(row=1,column=2)

txtShakes=Entry(f1,font=('arial',16,'bold'),textvariable=shakes,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtShakes.grid(row=1,column=3)

lblcostmeal=Label(f1,font=('arial',16,'bold'),text='Cost of Meal',bd=16,anchor='w')
lblcostmeal.grid(row=2,column=2)

txtcostmeal=Entry(f1,font=('arial',16,'bold'),textvariable=meal,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtcostmeal.grid(row=2,column=3)

lblservicecharge=Label(f1,font=('arial',16,'bold'),text='GST',bd=16,anchor='w')
lblservicecharge.grid(row=3,column=2)

txtservicecharge=Entry(f1,font=('arial',16,'bold'),textvariable=gst,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtservicecharge.grid(row=3,column=3)

lblsubtotal=Label(f1,font=('arial',16,'bold'),text='Sub Total',bd=16,anchor='w')
lblsubtotal.grid(row=4,column=2)

txtsubtotal=Entry(f1,font=('arial',16,'bold'),textvariable=subtotal,bd=10,bg="powder blue",justify="right",insertwidth=4)
txtsubtotal.grid(row=4,column=3)

lbltotal=Label(f1,font=('arial',16,'bold'),text='Total',bd=16,anchor='w')
lbltotal.grid(row=5,column=2)

txttotal=Entry(f1,font=('arial',16,'bold'),textvariable=total,bd=10,bg="powder blue",justify="right",insertwidth=4)
txttotal.grid(row=5,column=3)
'''================================================Button=========================================='''
btntotal=Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',16,"bold"),width=10,text="Total",bg="powder blue",command=ref)
btntotal.grid(row=7,column=1)

btnreset=Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',16,"bold"),width=10,text="Reset",bg="powder blue",command=reset)
btnreset.grid(row=7,column=2)

btnexit=Button(f1,padx=16,pady=8,bd=16,fg="red2",font=('arial',16,"bold"),width=10,text="Exit",bg="powder blue",command=qexit)
btnexit.grid(row=7,column=3)
window.mainloop()

import tkinter

root = tkinter.Tk()
values = "0"
inputs = tkinter.StringVar()
inputs.set(values)      # 초기값 설정한 0 보여줌


# 입력받는 함수
def inputValue(any):
    global values
    values = values + str(any)
    if(values[0]=='0'):
        if (values[1]=='+' or values[1]=='-' or values[1]=='*' or values[1]=='/' or values[1]=='.'):
            values = values[0:]
        else:
            values = values[1]
    else:
        pass
    inputs.set(values)

# 사칙연산 기호(+,-,*,/) 처리 함수(연속 입력 안되게...)
def inputSign(any):
    global values
    if (values[-1]=='+' or values[-1]=='-' or values[-1]=='*' or values[-1]=='/'):
        values = values[:-1] + str(any)
        inputs.set(values)
    elif(values[-1]=='.'):
        values = values[:-1] + str(any)
        inputs.set(values)
    else:
        values = values + str(any)
        inputs.set(values)




# 잘못 입력시 제거하는 함수 (←)
def delete(any):
    global values
    values = values[:-1]
    inputs.set(values)

# 입력 값 초기화(C)
def clear(any):
    global values
    values="0"
    inputs.set(values)

# 계산
def calculate(any):
    try:
        global values
        result = str(eval(values))
        inputs.set(result)
        values=result
    except:
        inputs.set("Error")
        values=""


# GUI 구현
# textvariable - Entry 위젯의 문자열 저장하는 Tk()변수(StringVar)
root.geometry("205x220")    # size 지정
root.resizable(0,0)         # 창 고정
display = tkinter.Entry(root, textvariable = inputs, font=('arial',15,'bold'), justify='right',width=15)
display.grid(row=1,column=0, columnspan=4,ipadx=15,ipady=10)
# shows = tkinter.Entry(display, textvariable = inputs,state = 'readonly',font = ('arial',10), justify='right',width=25)
# shows.grid(row=0,column=0, columnspan=4)

btnBack=tkinter.Button(root,text='BackSpace',command=lambda:delete('←'),width=5,height=2)
btnBack.grid(row=2, column=0,columnspan=2,ipadx=25)
btnCls=tkinter.Button(root,text='Clear',command=lambda:clear('C'),width=5,height=2)
btnCls.grid(row=2, column=2,columnspan=2,ipadx=25)

btn7=tkinter.Button(root,text='7',command=lambda:inputValue('7'),width=5,height=2)
btn7.grid(row=3, column=0)
btn8=tkinter.Button(root,text='8',command=lambda:inputValue('8'),width=5,height=2)
btn8.grid(row=3, column=1)
btn9=tkinter.Button(root,text='9',command=lambda:inputValue('9'),width=5,height=2)
btn9.grid(row=3, column=2)
btnDiv=tkinter.Button(root,text='/',command=lambda:inputSign('/'),width=5,height=2)
btnDiv.grid(row=3, column=3)

btn4=tkinter.Button(root,text='4',command=lambda:inputValue('4'),width=5,height=2)
btn4.grid(row=4, column=0)
btn5=tkinter.Button(root,text='5',command=lambda:inputValue('5'),width=5,height=2)
btn5.grid(row=4, column=1)
btn6=tkinter.Button(root,text='6',command=lambda:inputValue('6'),width=5,height=2)
btn6.grid(row=4, column=2)
btnMul=tkinter.Button(root,text='*',command=lambda:inputSign('*'),width=5,height=2)
btnMul.grid(row=4, column=3)

btn1=tkinter.Button(root,text='1',command=lambda:inputValue('1'),width=5,height=2)
btn1.grid(row=5, column=0)
btn2=tkinter.Button(root,text='2',command=lambda:inputValue('2'),width=5,height=2)
btn2.grid(row=5, column=1)
btn3=tkinter.Button(root,text='3',command=lambda:inputValue('3'),width=5,height=2)
btn3.grid(row=5, column=2)
btnMinus=tkinter.Button(root,text='-',command=lambda:inputSign('-'),width=5,height=2)
btnMinus.grid(row=5, column=3)

btn0=tkinter.Button(root,text='0',command=lambda:inputValue('0'),width=5,height=2)
btn0.grid(row=6, column=0)
btnPoint=tkinter.Button(root,text='.',command=lambda:inputSign('.'),width=5,height=2)
btnPoint.grid(row=6, column=1)
btnPlus=tkinter.Button(root,text='+',command=lambda:inputSign('+'),width=5,height=2)
btnPlus.grid(row=6, column=2)
btnEqual=tkinter.Button(root,text='=',command=lambda:calculate('='),width=5,height=2)
btnEqual.grid(row=6, column=3)


root.mainloop()

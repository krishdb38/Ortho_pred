import tkinter as tk

def evaluate(event):
    res.configure(text = "Result: " + str(eval(entry.get())))

w = tk.Tk()
tk.Label(w, text="Your Expression:").pack()
entry = tk.Entry(w)
entry.bind("", evaluate)
entry.pack()
res = tk.Label(w)
res.pack()
w.mainloop()
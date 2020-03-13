from alive_progress import alive_bar
items  = range(100000000000)
with alive_bar(len(items)) as bar:
    for item in items:
        i = 0
        i +=2
        bar()
import bs4
import requests
def getQuote():
    try:
        res = requests.get("https://www.brainyquote.com/quotes_of_the_day.html")
        soup = bs4.BeautifulSoup(res.text, 'lxml')
        quote = soup.find('img', {"class":"p-qotd"})

        text = quote['alt']
        #print(text)

        return text
    except ConnectionError as e:
        print("Connection error, please check internet connection", e)
print(getQuote())
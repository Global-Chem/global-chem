
from bs4 import BeautifulSoup

if __name__ == '__main__':

    soup = BeautifulSoup(open('test.html', 'r'), features="lxml")

    tables = soup.findAll('table', attrs={'class':'finder'})
    total_list = []

    for table in tables:

        elements = table.findAll('td')

        for i in elements:

            text = ''.join(i.renderContents().decode('utf-8').split()[1:])
            if len(text) != 0:
                total_list.append(text)

    for element in total_list:

        element = element.replace('<em>', '').replace(r'</em>', '').lower()
        print (f"'{element}': '',")




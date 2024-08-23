import requests
from bs4 import BeautifulSoup
import os

def get_links_from_table(url):
    # Fetch the HTML content
    response = requests.get(url)
    html_content = response.text

    # Parse the HTML
    soup = BeautifulSoup(html_content, 'html.parser')

    # Find the first table in the body
    table = soup.body.find('table')

    # Extract links from the second column of each row
    links = []
    if table:
        rows = table.find_all('tr')
        for row in rows:
            columns = row.find_all('td')
            if len(columns) >= 2:
                second_column = columns[1]
                anchors = second_column.find_all('a', href=True)
                for anchor in anchors:
                    links.append(anchor['href'])

    return links

url_root = "https://www.globsnow.info/swe/archive_v3.0/L3A_daily_SWE/NetCDF4/"
result = get_links_from_table(url_root)
for yy in result:
    os.system("wget " + url_root + yy + " -O " + os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'results', 'GlobSnow', yy))
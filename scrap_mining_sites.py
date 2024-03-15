import requests
from bs4 import BeautifulSoup, NavigableString
import re
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from urllib.parse import urljoin
import pandas as pd
import numpy as np
from tqdm import tqdm
import os


def dms_to_decimal(dms_str):
    # Regular expression to match DMS components in the string
    dms_str = re.sub(r'[a-zA-Z\s()\.]', '', dms_str)

    degree, rest = dms_str.split("°")
    degree = float(degree)
    rest = rest.split("'")
    if len(rest[0]) > 1:
        minute = float(rest[0])
    else:
        minute = 0
    if len(rest) > 1:
        rest = rest[1].replace("'", "")
        second = float(rest)
    else:
        second = 0
    result = degree + minute / 60 + second / 3600
    return result


options = Options()
options.headless = True
driver = webdriver.Chrome(executable_path = "C:/Users/ywo/OneDrive - Oak Ridge National Laboratory/Projects/2023 ELM-ERW LDRD/data/chromedriver-win64/chromedriver.exe", options=options)


def get_page(aurl):
    driver.get(aurl)
    html = driver.page_source
    return html


main_url = 'https://www.mindat.org/locentries.php?p=3366&m=48492'
html = get_page(main_url)
soup = BeautifulSoup(html, 'html.parser')


# Regular expression to match 'loc-<bunch of numbers>.html'
pattern = re.compile('loc-\d+\.html')


# Find all <a> tags with href matching the pattern
matching_links = [a['href'] for a in soup.find_all('a', href=pattern)]
site_names = [a.get_text() for a in soup.find_all('a', href=pattern)]


# Print the matching links
site_list = pd.DataFrame(np.nan, index = site_names, columns = ['lat', 'lon'])
for site, link in tqdm(zip(site_names, matching_links)):
    site_url = urljoin("https://www.mindat.org/", link)
    site_html = get_page(site_url)
    soup = BeautifulSoup(site_html, 'html.parser')
    google_maps_link = soup.find('a', href = re.compile('.*maps\.google\.com.*'))
    if google_maps_link:
        lat, lon = google_maps_link.get_text().split(",")
        site_list.loc[site, 'lat'] = float(lat)
        site_list.loc[site, 'lon'] = float(lon)
    else: 
        # try falling back to using the WGS 84 formatted coords
        divs = soup.find_all('div', class_ = "LFtr")

        # Iterate through the <div> elements
        for i in range(len(divs) - 1):
            # Check if the content of the current div is "Latitude & Longitude (WGS84):"
            if "Latitude & Longitude (WGS84):" in divs[i].get_text():
                # The next div is the one we are interested in
                target_div = divs[i].find_all('div', class_ = "LFtd")[0]
                wgs84 = ''
                for content in target_div:
                    if isinstance(content, NavigableString):
                        wgs84 += content.strip()

                north, west = wgs84.split(",")
                lat = dms_to_decimal(north)
                lon = - dms_to_decimal(west)
                site_list.loc[site, 'lat'] = float(lat)
                site_list.loc[site, 'lon'] = float(lon)

                break

        # if WGS 84 doesn't exist, this is not a site

site_list.to_csv(os.path.join('..', 'data', 'mining_sites.csv'))

site_list.dropna(axis = 0, how = 'all').to_csv(os.path.join('..', 'data', 'mining_sites_valid.csv'))

driver.quit()
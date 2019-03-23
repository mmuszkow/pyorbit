#!/usr/bin/env python3

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from time import sleep
from datetime import datetime
import sys

def dump():
    norad_ids = [20253, 17181] #[40069, 33591, 25338, 28654]
    #driver = webdriver.Safari()
    #driver = webdriver.PhantomJS()
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    driver = {}
    for norad_id in norad_ids:
        driver[norad_id] = webdriver.Chrome(options=chrome_options)
        driver[norad_id].get('http://www.n2yo.com/?s=%d' % norad_id)
    while True:
        for norad_id in norad_ids:
            tries = 0
            try:
                lon = lat = utc = alt = None
                while not lon or not lat or not utc or not alt:
                    if tries > 60: raise
                    tries += 1
                    sleep(1)
                    utc = driver[norad_id].find_element_by_id(id_='utctime').text
                    if utc: utc = datetime.utcnow().strftime('%Y-%m-%d_') + utc
                    lon = driver[norad_id].find_element_by_id(id_='satlng').text
                    lat = driver[norad_id].find_element_by_id(id_='satlat').text
                    alt = driver[norad_id].find_element_by_id(id_='sataltkm').text
                print(norad_id, utc, lon, lat, alt, flush=True)
            except Exception as e:
                print(e, file=sys.stderr)
                try: driver[norad_id].close()
                except: pass
                try: driver[norad_id] = webdriver.Chrome(chrome_options=chrome_options)
                except Exception as e: print(e, file=sys.stderr)
                sleep(1)
        
        sleep(1)

dump()


'''
OMIM api don't have moi. Scrape the website for it
'''
import json
from collections import Counter
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.firefox.options import Options
from selenium.common.exceptions import NoSuchElementException, StaleElementReferenceException
from lib import helper

remote = True
baseurl = 'https://www.omim.org/entry/'


def main(params):
    with open(params['simulation']['annotated_variants'], 'rt') as inf, open(params['simulation']['omim_moi'], 'wt') as outf:
        data = json.load(inf)
        omims = set()
        for gene in data:
            omim = Counter(gene['clinvar']['omim']).most_common(1)[0][0]
            omims.add(omim)

        if remote:
            from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
            driver = webdriver.Remote(
                command_executor='http://127.0.0.1:4444/wd/hub', desired_capabilities=DesiredCapabilities.FIREFOX)
        else:
            driver = webdriver.Firefox()
        wait = WebDriverWait(driver, 10)

        for omim in omims:
            print(omim)
            url = baseurl + omim
            driver.get(url)
            # dismiss donate
            #donate = driver.find_element_by_id('donationPopupModal')
            # if donate:
            #    donate.find_element_by_id('donationPopupCancel').click()
            moi_index = None
            gene_index = None
            try:
                table = driver.find_element_by_tag_name('table')
                header = table.find_element_by_tag_name(
                    'thead').find_elements_by_tag_name('th')
                entries = table.find_element_by_tag_name(
                    'tbody').find_elements_by_tag_name('tr')
            except NoSuchElementException:
                print('no gene or moi found for this omim')
                continue

            # get moi_index and gene_index, relative to len(header)
            # since there might be a td in the front spanning multiple rows
            header = [i.get_attribute('innerHTML').strip() for i in header]
            print(header)
            try:
                moi_index = header.index('Inheritance') - len(header)
            except ValueError:
                print('no visible table, pass')
                continue
            try:
                gene_index = header.index('Gene/Locus') - len(header)
            except ValueError:
                pass
            for tr in entries:
                moi = ''
                gene = ''
                try:
                    moi = tr.find_elements_by_tag_name('td')[moi_index].find_element_by_tag_name(
                        'abbr').get_attribute('innerHTML').strip()
                except NoSuchElementException:
                    pass
                if gene_index is not None:
                    try:
                        gene = tr.find_elements_by_tag_name('td')[gene_index].find_element_by_tag_name(
                            'span').get_attribute('innerHTML').strip()
                    except NoSuchElementException:
                        pass
                outf.write('\t'.join([omim, moi, gene]) + '\n')


if __name__ == '__main__':
    config_file = 'configure.cfg'
    config = helper.parse_config(config_file)
    main(config)
    print('===done===')

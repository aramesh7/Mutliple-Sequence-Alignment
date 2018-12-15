import unittest
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
import pandas
import msa_utils as util
import os


class TestAlignmentFunctions(unittest.TestCase):
    def __init__(self):
        super(TestAlignmentFunctions, self).__init__()
        test_url = "http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh"
        chrome_options = Options()
        prefs = {"download.default_directory": "/Users/Ashwin/Desktop/CS466/Multiple-Sequence-Alignment"}
        chrome_options.add_experimental_option("prefs", prefs)
        # chrome_options.add_argument("--headless")
        self.driver = webdriver.Chrome(options=chrome_options)
        self.driver.get(test_url)

    def test_sequence_sequence(self):
        data = ["AC", "ACTG\n"]

        for i, e in enumerate(['sequence_1', 'sequence_2', ]):
            try:
                element = WebDriverWait(self.driver, 5).until(EC.presence_of_element_located((By.ID, e)))
                element.clear()
                element.send_keys(data[i])
            except TimeoutException:
                print("Loading took too much time!")

        for i, e in enumerate(["table_vertical_download", "table_download", "table_horizontal_download"]):
            try:
                element = WebDriverWait(self.driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, e)))
                element.click()
            except TimeoutException:
                print("Downloading took too much time!")

        paths = WebDriverWait(self.driver, 120, 1).until(every_downloads_chrome)
        for i, e in enumerate(["table.csv", "table (1).csv", "table (2).csv"]):
            os.rename(e, file_names[i] + ".csv")

    def test_sequence_profile(self):
        pass

    def test_profile_profile(self):
        pass

    def every_downloads_chrome(self, driver):
        if not driver.current_url.startswith("chrome://downloads"):
            driver.get("chrome://downloads/")
        return driver.execute_script("""
            var items = downloads.Manager.get().items_;
            if (items.every(e => e.state === "COMPLETE"))
                return items.map(e => e.file_url);
            """)

if __name__ == '__main__':
    unittest.main()





def every_downloads_chrome(driver):
    if not driver.current_url.startswith("chrome://downloads"):
        driver.get("chrome://downloads/")
    return driver.execute_script("""
        var items = downloads.Manager.get().items_;
        if (items.every(e => e.state === "COMPLETE"))
            return items.map(e => e.file_url);
        """)

for i,e in enumerate(["table_vertical_download", "table_download", "table_horizontal_download"]):
    try:
        element = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, e)))
        element.click()
    except TimeoutException:
        print("Downloading took too much time!")

paths = WebDriverWait(self.driver, 120, 1).until(every_downloads_chrome)
for i,e in enumerate(["table.csv", "table (1).csv", "table (2).csv"]):
    os.rename(e, file_names[i]+".csv")

# MUSCLE, CLUSTALW
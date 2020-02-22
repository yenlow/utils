#Create a gspread class and extract the data from the sheets
#requires:
# 1. Google API credentials json_key file path
# 2. scope e.g. ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
# 3. gspread_url e.g. 'https://docs.google.com/spreadsheets/d/1itaohdPiAeniCXNlntNztZ_oRvjh0HsGuJXUJWET008/edit?usp=sharing'

import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas as pd

class gspread_obj(object):
    """
    Create a google spreadsheet instance to download sheet(s) and merge them
    Requires spreadsheet url and Google API json key file

    Examples:
        >>>> gc = gspread_obj()

        >>>> gc.login('home/user/google_api_key.json')

        >>>> gc.get_sheets('https://docs.google.com/spreadsheets/d/1itaohdPiAeniCXNlntNztZ_oRvjh0HsGuJXUJWET008/edit?usp=sharing')

        >>>> df = gc.merge_sheets()

    """
    def __init__(self):
        self.scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
        self.client = None # gspread.Client object
        self.sheets = None

    def login(self, credentials_google: str):
        #set Google spreadsheet credentials
        credentials = ServiceAccountCredentials.from_json_keyfile_name(credentials_google, self.scope)
        self.client = gspread.authorize(credentials)

    def get_sheets(self, gspread_url: str):
        #Get Google sheet instance
        wks = self.client.open_by_url(gspread_url)
        self.sheets = wks.worksheets()

    def merge_sheets(self):
        if self.sheets is None:
            print('No sheets are found!')
            df = None

        elif len(self.sheets)==1:
            data = self.sheets[0].get_all_values()
            header = data.pop(0)
            df = pd.DataFrame(data, columns=header)

        elif len(self.sheets)>1:
            #read all the sheets
            df_list = []
            for s in self.sheets:
                data = s.get_all_values()
                header = data.pop(0)
                df = pd.DataFrame(data, columns=header)
                df_list.append(df)
            df = pd.concat(df_list, axis=0, join='outer', sort=False)

        else:
            print("self.sheets must be a list of sheet(s)!")
            df = None

        if df is not None:
            print("Columns: ", df.columns)
            print("{} Rows x {} Columns".format(df.shape[0],df.shape[1]))
        return df

# Config params
import os

#for Google API
gspread_url = 'https://docs.google.com/spreadsheets/d/1itaohdPiAeniCXNlntNztZ_oRvjh0HsGuJXUJWET008/edit?usp=sharing'

def google_api():
    credentials_google = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS") #or the key string
    return credentials_google

#for Github v3 API
def github_api():
    credentials_github = os.environ.get("GITHUB_ACCESS_TOKEN") #or the key string
    return credentials_github

#for Twitter
def twitter_api():
    consumer_key = os.environ.get("TWITTER_CONSUMER_KEY")
    consumer_secret = os.environ.get("TWITTER_CONSUMER_SECRET")
    access_token = os.environ.get("TWITTER_ACCESS_TOKEN")
    access_secret = os.environ.get("TWITTER_ACCESS_SECRET")
    return consumer_key, consumer_secret, access_token, access_secret



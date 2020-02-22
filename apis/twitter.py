from apis.config import twitter_api
import pickle
import tweepy as tw
import pandas as pd

pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 2000)
pd.set_option('display.width', 1000)

#https://bhaskarvk.github.io/2015/01/how-to-use-twitters-search-rest-api-most-effectively./
consumer_key, consumer_secret, access_token, access_secret = twitter_api()
auth = tw.OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_secret)
api = tw.API(auth, wait_on_rate_limit=True, wait_on_rate_limit_notify=True)

# Define the search term and the date_since date as variables
#https://developer.twitter.com/en/docs/tweets/search/api-reference/get-search-tweets
searchQuery = "covid -filter:retweets"
date_since = "2020-02-20"
maxTweets = 300 # Some arbitrary large number
tweetsPerQry = 100  # this is the max the API permits
sinceId = None #Returns results after this ID (i.e. more recent)
max_id = -float('inf') #Returns results before this ID (i.e. older)

# filenames
f_tweets = 'tweets.pkl'
f_df_tweets = 'df_tweets.pkl'

print("Start downloading up to {0} tweets...".format(maxTweets))
tweetCount = 0
while tweetCount < maxTweets:
    try:
        if (max_id <= 0):
            if (not sinceId):
                new_tweets = api.search(q=searchQuery, count=tweetsPerQry, lang='en')
            else:
                new_tweets = api.search(q=searchQuery, count=tweetsPerQry, lang='en',
                                        since_id=sinceId)
        else:
            if (not sinceId):
                new_tweets = api.search(q=searchQuery, count=tweetsPerQry, lang='en',
                                        max_id=str(max_id - 1))
            else:
                new_tweets = api.search(q=searchQuery, count=tweetsPerQry, lang='en',
                                        max_id=str(max_id - 1),
                                        since_id=sinceId)
        if not new_tweets:
            print("No more tweets found")
            break

        #Save/dump everything ASAP
        with open(f_tweets, 'a+b') as f:
            pickle.dump(new_tweets, f)
        # for tweet in new_tweets:
        #     f.write(jsonpickle.encode(tweet._json, unpicklable=False) + '\n')

        tweetCount += len(new_tweets)
        max_id = new_tweets[-1].id
        print("Downloaded {0} tweets".format(tweetCount))

    except tw.TweepError as e:
        print("some error : " + str(e)) # Just exit if any error
        break

print("Downloaded {0} tweets total - saved to {1}".format(tweetCount, f_tweets))


# Load whole pickle tweetstream
tweets = []
with open(f_tweets, 'rb') as fr:
    try:
        while True:
            tweets += pickle.load(fr)
    except EOFError:
        pass


# Get the selected fields
tweets_sel = []
for tweet in tweets:
    tweets_sel.append([
                        tweet.id
                        ,tweet.text
                        ,tweet.user.id
                        ,tweet.user.name
                        ,tweet.user.screen_name
                        ,tweet.user.location
                        ,tweet.user.statuses_count
                        ,tweet.user.verified
                        ,tweet.favorite_count
                        ,tweet.favorited
                        ,tweet.retweet_count
                        ,tweet.entities
                        #,tweet.entities['hashtags']
                        #,tweet.entities['user_mentions']['screen_name']
                    ])

colnames = [
            'id'
            ,'text'
            ,'user_id'
            ,'user_name'
            ,'user_screen_name'
            ,'user_location'
            ,'user_statuses_count'
            ,'user_verified'
            ,'favorite_count'
            ,'favorited'
            ,'retweet_count'
            ,'entities'
            ]


df_tweets = pd.DataFrame(data=tweets_sel, columns=colnames)

with open(f_df_tweets, 'wb') as f:
    pickle.dump(df_tweets, f)

# dir(tweet)
# ['__class__',
#  '__delattr__',
#  '__dict__',
#  '__dir__',
#  '__doc__',
#  '__eq__',
#  '__format__',
#  '__ge__',
#  '__getattribute__',
#  '__getstate__',
#  '__gt__',
#  '__hash__',
#  '__init__',
#  '__init_subclass__',
#  '__le__',
#  '__lt__',
#  '__module__',
#  '__ne__',
#  '__new__',
#  '__reduce__',
#  '__reduce_ex__',
#  '__repr__',
#  '__setattr__',
#  '__sizeof__',
#  '__str__',
#  '__subclasshook__',
#  '__weakref__',
#  '_api',
#  '_json',
#  'author',
#  'contributors',
#  'coordinates',
#  'created_at',
#  'destroy',
#  'entities',
#  'favorite',
#  'favorite_count',
#  'favorited',
#  'geo',
#  'id',
#  'id_str',
#  'in_reply_to_screen_name',
#  'in_reply_to_status_id',
#  'in_reply_to_status_id_str',
#  'in_reply_to_user_id',
#  'in_reply_to_user_id_str',
#  'is_quote_status',
#  'lang',
#  'metadata',
#  'parse',
#  'parse_list',
#  'place',
#  'retweet',
#  'retweet_count',
#  'retweeted',
#  'retweeted_status',
#  'retweets',
#  'source',
#  'source_url',
#  'text',
#  'truncated',
#  'user']
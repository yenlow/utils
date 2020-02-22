# Download useful scripts and functions from github
# Most are in R
from apis.github import gh_obj
from apis.config import credentials_github

gh = gh_obj()
gh.login(credentials_github)

#from https://github.com/jameshay218/case_to_infection
files = ['code/date_functions.R',
         'code/analysis_functions.R',
         'code/augmentation_functions.R',
         'code/pull_and_clean_linelist.R',
         'code/analysis_standalone.R']

for f in files:
    gh.get_file(repo="jameshay218/case_to_infection", file=f, outdir="other_authors")

#from https://github.com/beoutbreakprepared/nCoV2019
files = ['entry-checker.R',
 #        'de_dupe_functions.R',
 #        'de_dupe_linelist.R',
         'code/processing_functions/routine_summary.R']

for f in files[-1]:
    gh.get_file(repo="beoutbreakprepared/nCoV2019", file=f, outdir="other_authors")
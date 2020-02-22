# 1. Download  scripts from public github (use wget or via Github V3 API)

from apis.github import Github


def get_rawfile_from_github(raw_url: str, outdir: str = None) -> None:
    """
    Download raw script from github. URL should point to the raw script off github

    Args:
        raw_url (str): url of raw script off github (should be publicly accessible)
        outdir (str): directory for saving the script; defaults to current path otherwise
    Example:
        >>> get_rawfile_from_github('https://raw.githubusercontent.com/beoutbreakprepared/nCoV2019/master/entry-checker.R','other_authors')

    """
    from wget import download, filename_from_url
    filename = filename_from_url(raw_url)
    if outdir:
        outpath = outdir + '/' + filename
    else:
        outpath = filename
    download(raw_url,outpath)


class gh_obj(object):
    """
    Creates a github instance to download file from repo via the Github V3 API
    This requires Github V3 API access token
    See https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line

    Examples:
        >>>gh = gh_obj("jameshay218/case_to_infection")

        >>>gh.login(<accesstoken>)

        >>>gh.get_file(file="code/analysis_functions.R", outdir="other_authors")
    """
    def __init__(self):
        self.instance = None

    def login(self, credentials_github: str):
        self.instance = Github(credentials_github)

    def get_file(self, repo: str, file: str, outdir: str = None) -> None:
        if not repo:
            raise ValueError('Please provide repo as str')

        repoinstance = self.instance.get_repo(repo)
        contents = repoinstance.get_contents(file)

        if outdir:
            outpath = outdir + "/" + contents.name
        else:
            outpath = contents.name

        f = open(outpath, "wb")
        f.write(contents.decoded_content)
        f.close()

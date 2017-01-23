import urllib2


class Download:
    def __init__(self, url='', verbose=False):
        self.url = url
        self.verbose = verbose
        self.html = ''
        self.error = False

    def download(self, user_agent='wswp', num_retries=2):
        self.html, self.error = download(self.url, self.verbose,
                                         user_agent=user_agent,
                                         num_retries=num_retries)


def download(url, verbose, user_agent='wswp', num_retries=2):
    """
    Function to download contents from a given url

    Input:
            url: str
            string with the url to download from

            user_agent: str
            Default 'wswp'

            num_retries: int
            Number of times to retry downloading
            if there is an error

            verbose: bool
            Print out url and errors

    Output:
            returns: str
            string with contents of given url
    """

    error = False
    if verbose:
        print('Downloading:', url)
    headers = {'User-agent': user_agent}
    request = urllib2.Request(url, headers=headers)
    try:
        html = urllib2.urlopen(request).read()
    except urllib2.URLError as e:
        if verbose:
            print('Download error:', e.reason)
        html = None
        if num_retries > 0:
            if hasattr(e, 'code') and 500 <= e.code < 600:
                # retry 5XX HTTP errors
                return download(url, user_agent, num_retries - 1)[0]
            # elif hasattr(e, 'code') and e.code == 404:
            else:
                error = True

    return html, error

from urllib import request, error


class Download:
    def __init__(self, url='', verbose=False, timeout=5):
        self.url = url
        self.verbose = verbose
        self.html = ''
        self.error = False
        self.timeout = timeout

    def download(self, user_agent='wswp', num_retries=2):
        # self.html, self.error = download(self.url, self.verbose, user_agent=user_agent, num_retries=num_retries,
        #                                  timeout=self.timeout)
        try:
            self.html = download(self.url, self.verbose, user_agent=user_agent, num_retries=num_retries,
                                 timeout=self.timeout)
        except IOError:
            raise ValueError("Could not download url.")


def download(url, verbose, user_agent='wswp', num_retries=2, decoding_format='utf-8', timeout=5):
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

            decoding: "utf-8"

    Output:
            returns: str
            string with contents of given url
    """

    #html_error = False
    if verbose:
        print('Downloading:', url)
    headers = {'User-agent': user_agent}
    request_obj = request.Request(url, headers=headers)
    try:
        with request.urlopen(request_obj, timeout=timeout) as response:
            html = response.read()
    except error.URLError as e:
        if verbose:
            print('Download error:', e.reason)
        # html = None
        # if num_retries > 0:
        #     if hasattr(e, 'code') and 500 <= e.code < 600:
        #         # retry 5XX HTTP errors
        #         return download(url, user_agent, num_retries - 1)[0]
        #     # elif hasattr(e, 'code') and e.code == 404:
        #     else:
        #         html_error = True
        raise IOError(e.reason)

    return html.decode(decoding_format)

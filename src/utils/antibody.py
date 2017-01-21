#from downloads import Download
import re
import urllib2

class Antibody:
    """
    TODO: write description
    """

    def __init__(self, sequence='', name='', numbering=None):
        self.raw_sequence = sequence
        self.sequence = sequence.replace('-','')
        self.name = name
        self.numbering = numbering
        self.url=''

    def apply_numbering(self, server='abysis', numbering_scheme='chothia'):

        assert numbering_scheme.lower() in ['chothia', 'chothia_ext', 'kabat']

        # check which server to use to get numbering
        if server.lower() == 'abysis':
            # find out which numbering scheme to use
            if numbering_scheme.lower() == 'chothia':
                scheme = '-c'
            elif numbering_scheme.lower() == 'chotia_ext':
                scheme = '-a'
            else:
                scheme = '-k'

            self.url = 'http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq={}&scheme={}'.format(self.sequence,
                                                                                                      scheme)

            numbering_table, error = download(self.url, verbose=False)

            parsed_numbering_table = re.findall('[\S| ]+', numbering_table)

            self.numbering = [x[:-2] for x in parsed_numbering_table]

            # TODO: add more server options


def download(url, user_agent='wswp', num_retries=2, verbose=True):
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
        print 'Downloading:', url
    headers = {'User-agent': user_agent}
    request = urllib2.Request(url, headers=headers)
    try:
        html = urllib2.urlopen(request).read()
    except urllib2.URLError as e:
        if verbose:
            print 'Download error:', e.reason
        html = None
        if num_retries > 0:
            if hasattr(e, 'code') and 500 <= e.code < 600:
                # retry 5XX HTTP errors
                return download(url, user_agent, num_retries - 1)[0]
            # elif hasattr(e, 'code') and e.code == 404:
            else:
                error = True

    return html, error
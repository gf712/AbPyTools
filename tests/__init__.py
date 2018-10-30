from urllib import request

ABNUM_URL = 'http://www.bioinf.org.uk/abs/abnum'
IGBLAST_URL = 'https://www.ncbi.nlm.nih.gov/igblast/'


# Helper functions
def check_connection(URL, timeout=5):
    try:
        request.urlopen(url=URL, timeout=timeout)
        return True
    except request.URLError:
        return False


def read_sequence(path):
    with open(path, 'r') as f:
        data = f.readlines()[1]
    return data


def read_sequence_from_file(file):
    with open(file, 'r') as f:
        data = f.readlines()[0]
    return data

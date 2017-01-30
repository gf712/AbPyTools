import sys


class PythonConfig:
    def __init__(self):
        self.backend = None

    def get_ipython_info(self):
        self.backend = ipython_info()


def ipython_info():
    # code obtained from stackoverflow forum:
    # http://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    ip = 'other'
    if 'ipykernel' in sys.modules:
        ip = 'notebook'
    elif 'IPython' in sys.modules:
        ip = 'terminal'
    return ip

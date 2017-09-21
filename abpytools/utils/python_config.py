import sys
from matplotlib import pyplot as plt


class PythonConfig:
    def __init__(self):
        self._backend = get_ipython_info()
        self._matplotlib_interactive = plt.isinteractive()

    @property
    def ipython_info(self):
        return self._backend

    @property
    def matplotlib_interactive(self):
        return self._matplotlib_interactive


def get_ipython_info():
    # code obtained from stackoverflow forum:
    # http://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    ip = 'other'
    if 'ipykernel' in sys.modules:
        ip = 'notebook'
    elif 'IPython' in sys.modules:
        ip = 'terminal'
    return ip

import os


class Home():

    def __init__(self):

        self.homedir = os.path.dirname(os.path.abspath(__file__))

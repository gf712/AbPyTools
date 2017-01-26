import json
from home import Home

# TODO: sort out a single method of loading data
# - best would be to put everything in json files with appropriate names
# - or put everything in a single file -> could become confusing to edit

class DataLoader:

    def __init__(self, numbering='', amino_acid_property=[], misc=[]):

        self.numbering = numbering
        self.amino_acid_property = amino_acid_property
        self.misc = misc
        self.directory_name = Home().homedir

    def get_data(self):
        # type: () -> list (name of positions)
        if len(self.numbering) > 0:
            with open('{}/data/NumberingSchemes/{}.txt'.format(self.directory_name, self.numbering), 'r') as f:
                data = [x.replace('\n', '') for x in f.readlines()]

        elif len(self.misc) > 0:
            with open('{}/data/{}'.format(self.directory_name, self.misc[0]), 'r') as f:
                data = json.load(f)[self.misc[1]][self.misc[2]]

        elif len(self.amino_acid_property) > 0:
            measurement = self.amino_acid_property[0]
            database = self.amino_acid_property[1]
            with open('{}/data/AminoAcidProperties.json'.format(self.directory_name), 'r') as f:
                data = json.load(f)[measurement][database]

        else:
            data = ''

        return data

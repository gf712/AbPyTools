import json

from home import Home


class DataLoader:

    def __init__(self, numbering='', amino_acid_property=[]):

        self.numbering = numbering
        self.amino_acid_property = amino_acid_property
        self.directory_name = Home().homedir

    def get_data(self):
        # type: () -> list (name of positions)
        if len(self.numbering) > 0:
            with open('{}/data/NumberingSchemes/{}.txt'.format(self.directory_name, self.numbering), 'r') as f:
                data = [x.replace('\n', '') for x in f.readlines()]
        elif len(self.amino_acid_property) > 0:
            measurement = self.amino_acid_property[0]
            database = self.amino_acid_property[1]

            with open('{}/data/AminoAcidProperties.json'.format(self.directory_name), 'r') as f:
                data = json.load(f)[measurement][database]
        else:
            data = ''
        return data

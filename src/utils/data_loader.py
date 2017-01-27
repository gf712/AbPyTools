import json
from home import Home

# TODO: sort out a single method of loading data
# - best would be to put everything in json files with appropriate names
# - or put everything in a single file -> could become confusing to edit


class DataLoader:

    def __init__(self, data=[], data_type='', amino_acid_property=[]):

        """

        :param data: list with the information required to access data in json file
        :param data_type: str must be one of the following: 'CDR_positions', 'NumberingSchemes' or 'AminoAcidProperties'
        :param amino_acid_property: temporary parameter
        :param misc: temporary parameter
        """
        self.amino_acid_property = amino_acid_property
        self.directory_name = Home().homedir
        self.data = data
        self.data_types = ['CDR_positions', 'NumberingSchemes', 'AminoAcidProperties']

        if data_type not in self.data_types:
            raise ValueError("{} is not a valid data type. Available data types: {}".format(
                data_type, ' ,'.join(self.data_types)))
        else:
            self.data_type = data_type

        # checks values when object is instantiated
        if self.data_type == 'CDR_positions' or self.data_type == 'NumberingSchemes':
            if len(self.data) != 2:
                raise ValueError("Expected 2, instead of {} values.".format(len(self.data)))
            if self.data[0] not in ['chothia']:
                raise ValueError("Got {}, but only {} is available".format(self.data[0], 'chothia'))
            if self.data[1] not in ["light", "heavy"]:
                raise ValueError("Got {}, but only light and heavy are available".format(self.data[1]))
        else:
            if len(self.data) != 2:
                raise ValueError("Expected 2, instead of {} values.".format(len(self.data)))
            if self.data[0] not in ["hydrophobicity", "pI", "MolecularWeight"]:
                raise ValueError("Got {}, but only {} are available".format(self.data[0],
                                                                            "hydrophobicity, pI and MolecularWeight"))

            if self.data[1] not in ["kdHydrophobicity", "wwHydrophobicity", "hhHydrophobicity", "mfHydrophobicity",
                                    "ewHydrophobicity", "EMBOSS", "DTASetect", "Solomon", "Sillero", "Rodwell",
                                    "Wikipedia", "Lehninger", "Grimsley", "average", "monoisotopic"]:
                raise ValueError("Got {}, but only {} are available".format(self.data[1],
                                                                            """
                                                                            kdHydrophobicity ,wwHydrophobicity,
                                                                            hhHydrophobicity, mfHydrophobicity,
                                                                            ewHydrophobicity, EMBOSS, DTASetect,
                                                                            Solomon, Sillero, Rodwell,
                                                                            Wikipedia, Lehninger, Grimsley,
                                                                            average and monoisotopic
                                                                            """
                                                                            ))

    def get_data(self):

        if self.data_type == 'CDR_positions':
            with open('{}/data/CDR_positions.json'.format(self.directory_name), 'r') as f:
                # need to access numbering scheme and chain type
                data = json.load(f)[self.data[0]][self.data[1]]
        elif self.data_type == 'NumberingSchemes':
            with open('{}/data/NumberingSchemes.json'.format(self.directory_name), 'r') as f:
                # need to access numbering scheme and chain type
                data = json.load(f)[self.data[0]][self.data[1]]
        else:
            with open('{}/data/AminoAcidProperties.json'.format(self.directory_name), 'r') as f:
                data = json.load(f)[self.data[0]][self.data[1]]

        return data

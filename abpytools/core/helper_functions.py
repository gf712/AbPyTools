from ..utils import DataLoader
import itertools
import pandas as pd


available_regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']


def numbering_table_sequences(region, numbering_scheme, chain):
    # next two conditionals are used to only extract the needed data with FR and CDR positions
    # it makes the loading of the data quicker when there are only CDRs or FRs
    if any([True if x.startswith('CDR') else False for x in region]):
        cdr_list = DataLoader(data_type='CDR_positions',
                              data=[numbering_scheme, chain])
        cdrs = cdr_list.get_data()
    else:
        cdrs = {}

    if any([True if x.startswith('FR') else False for x in region]):
        fr_list = DataLoader(data_type='Framework_positions',
                             data=[numbering_scheme, chain])
        frs = fr_list.get_data()
    else:
        frs = {}

    # pack it all up into a single dictionary
    whole_sequence_dict = {**cdrs, **frs}
    # get the sequence list in the correct order (since region has been sorted before)
    whole_sequence_list = [whole_sequence_dict[x] for x in region]
    # unpack whole_sequence_list into a single list
    whole_sequence = list(itertools.chain.from_iterable(whole_sequence_list))

    return whole_sequence_dict, whole_sequence


def numbering_table_region(region):

    # if 'all' is chosen then region becomes a list with all
    if region == 'all':
        region = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']

    # if region is a string (i.e. 'CDR1') it becomes a list and in the next block of code is checked if it
    # is a valid selection
    if isinstance(region, str):
        region = [region]

    # checks if the selected regions are all in the available_regions
    if not set(region).issubset(set(available_regions)):
        raise ValueError("The chosen region is not available"
                         "Currently available regions: {}".format(available_regions))

    # make sure the regions are in a logical order
    region.sort(key=lambda x: available_regions.index(x))

    return region


def numbering_table_multiindex(region, whole_sequence_dict):
    # [[(CDR1, L23), (CDR1, L24),..], ..., [(CDR3, L100), ...]]
    # here we get the list of lists from the line above
    pre_region_map = [[(region_i, numbering) for numbering in whole_sequence_dict[region_i]] for region_i in
                      region]

    # which can be easily unpacked into a single list
    # region_map is a list of tuples that can be interpreted by pd.MultiIndex to form a two layer column system
    region_map = list(itertools.chain.from_iterable(pre_region_map))
    multi_index = pd.MultiIndex.from_tuples(tuples=region_map, names=['Region', 'Numbering'])

    return multi_index


def germline_identity_pd(heavy_identity, light_identity, names):

    h_germline_pd = pd.DataFrame(heavy_identity).T
    l_germline_pd = pd.DataFrame(light_identity).T

    l_columns = pd.MultiIndex.from_tuples([('Light', x) for x in l_germline_pd.columns], names=['Chain', 'Region'])
    h_columns = pd.MultiIndex.from_tuples([('Heavy', x) for x in h_germline_pd.columns], names=['Chain', 'Region'])
    average_columns = pd.MultiIndex.from_tuples([('Average', x) for x in l_germline_pd.columns],
                                                names=['Chain', 'Region'])

    # h = pd.DataFrame(index=h_germline_pd.index.tolist(),
    #                  columns=l_columns)
    # l = pd.DataFrame(index=l_germline_pd.index.tolist(),
    #                  columns=h_columns)
    #
    # l = l.apply(lambda x: l_germline_pd.loc[x.name], axis=1)
    # h = h.apply(lambda x: h_germline_pd.loc[x.name], axis=1)

    average = (h_germline_pd.as_matrix() + l_germline_pd.as_matrix()) / 2

    l_germline_pd.columns = l_columns
    h_germline_pd.columns = h_columns
    average = pd.DataFrame(average, columns=average_columns, index=names)

    l_germline_pd.index = names
    h_germline_pd.index = names

    return pd.concat([l_germline_pd, h_germline_pd, average], axis=1)

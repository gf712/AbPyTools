from ..utils import DataLoader
import itertools
import pandas as pd
import numpy as np


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


def germline_identity_pd(heavy_identity, light_identity, internal_heavy, internal_light, names):

    regions = ['CDR1', 'CDR2', 'CDR3', 'FR1', 'FR2', 'FR3', 'Total']

    columns = pd.MultiIndex.from_tuples([('Light', x) for x in regions] +
                                        [('Heavy', x) for x in regions] +
                                        [('Average', x) for x in regions],
                                        names=['Chain', 'Region'])

    df = pd.DataFrame(columns=columns, index=names)

    for column in columns:

        if column[0] == 'Light':
            df[column] = [light_identity[x][column[1]] if column[1] in light_identity[x] else np.NaN
                          for x in internal_light]
            # df[column] = list(map(lambda x: light_identity[x][column[1]] if column[1] in light_identity[x] else np.NaN,
            #                       internal_light))
        elif column[0] == 'Heavy':
            df[column] = [heavy_identity[x][column[1]] if column[1] in heavy_identity[x] else np.NaN
                          for x in internal_heavy]
            # df[column] = list(map(lambda x: heavy_identity[x][column[1]] if column[1] in heavy_identity[x] else np.NaN,
            #                       internal_heavy))
        else:
            df[column] = (df[('Light', column[1])] + df[('Heavy', column[1])]) / 2

    return df


def to_numbering_table(as_array, region, chain, heavy_chains_numbering_table,
                       light_chains_numbering_table, names, **kwargs):

    if chain == 'both':

        if as_array:
            t_heavy = heavy_chains_numbering_table(as_array=True, region=region, **kwargs)
            t_light = light_chains_numbering_table(as_array=True, region=region, **kwargs)

            data = np.concatenate((t_light, t_heavy), axis=1)

        else:
            t_heavy = heavy_chains_numbering_table(as_array=False, region=region, **kwargs)
            t_light = light_chains_numbering_table(as_array=False, region=region, **kwargs)

            t_heavy.reset_index(drop=True, inplace=True)
            t_light.reset_index(drop=True, inplace=True)

            data = pd.concat([t_light, t_heavy], axis=1, keys=['Light', 'Heavy'])

    elif chain == 'heavy':

        data = heavy_chains_numbering_table(as_array=as_array, region=region, **kwargs)

    elif chain == 'light':

        data = light_chains_numbering_table(as_array=as_array, region=region, **kwargs)

    else:
        raise ValueError("Unknown chain.")

    if not as_array:
        data.index = names

    return data

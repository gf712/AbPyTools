from abpytools import Chain
from ..utils import DataLoader

OPTIONS = ["CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3", "FR4"]


def expression_evaluator(expression, antibody_object):
    """
    A function to extract indices and amino position names information from a series of query commands
    Example: "index CDR3" returns index of amino acids in CDR3
             "amino acids CDR3" returns name of amino acids
             "numbering CDR3" returns the numbering of amino acids in CDR3 in the given sequence
             
    Params: expression
            antibody_object
    """

    if not (isinstance(antibody_object, Chain)):
        raise ValueError('Invalid antibody_object type')

    # first determine what to return
    if expression.startswith("index"):
        value = "index"
    elif expression.startswith("amino acid"):
        value = "amino acid"
    elif expression.startswith("numbering"):
        value = "numbering"
    else:
        raise ValueError("Invalid syntax")

    # remove the first part as it is not needed anymore
    expression = expression.replace(value, '').strip().split()

    # check if there is "exception"
    exception_expression = []
    if 'except' in expression:
        exception_expression_index = expression.index('except')
        exception_expression = expression[exception_expression_index+1:]
        expression = expression[:exception_expression_index]

    if 'all' in expression:
        expression = OPTIONS

    # evaluate remaining expression
    if set(expression).issubset(set(OPTIONS)) and set(exception_expression).issubset(set(OPTIONS)):
        # now extract the relevant values given what we know
        complete_numbering_dict = antibody_object.ab_regions()
        # need to unpack dictionaries into one (.ab_regions() returns a tuple (dict with CDRs, dict with FRs)
        complete_numbering_dict = {**complete_numbering_dict[0], **complete_numbering_dict[1]}
        index_query = [complete_numbering_dict[region] for region in expression]
        index_query_flat = [item for sublist in index_query for item in sublist]

        if value == 'index':
            return index_query_flat

        elif value == 'amino acid':
            return ''.join([antibody_object.sequence[i] for i in index_query_flat])

        else:
            data_loader = DataLoader(data_type='CDR_positions',
                                     data=[antibody_object.numbering_scheme, antibody_object.chain])
            whole_sequence_dict = data_loader.get_data()

            whole_sequence = [whole_sequence_dict[region] for region in expression]
            whole_sequence_flat = [item for sublist in whole_sequence for item in sublist]

            return [numbering_i for numbering_i in antibody_object.numbering if numbering_i in whole_sequence_flat]

    else:
        raise ValueError("""Invalid syntax, make sure you are selecting a valid region""")

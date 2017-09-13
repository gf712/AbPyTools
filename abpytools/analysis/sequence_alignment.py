from .analysis_helper_functions import load_alignment_algorithm, load_substitution_matrix


class SequenceAlignment:

    """
    Sequence alignment with two chain objects   
    """

    def __init__(self, target, collection, algorithm, substitution_matrix):

        """
        
        :param chain_1: 
        :param collection: 
        :param algorithm: 
        :param substitution_matrix: 
        """

        self._algorithm = algorithm
        self._substitution_matrix = load_substitution_matrix(substitution_matrix)
        self.target = target
        self._collection = collection
        self._aligned_collection = dict()
        self._alignment_scores = dict()

    def align_sequences(self, **kwargs):

        # perform the alignment for each chain object in collection and store results in dictionaries with keys
        # corresponding to names of the sequence to be aligned
        for seq in self._collection.antibody_objects:
            result = self._align(self.target, seq, **kwargs)
            self._aligned_collection[seq.name] = result[0]
            self._alignment_scores[seq.name] = result[1]

    def _align(self, seq_1, seq_2, **kwargs):

        # loads the function object that is then called
        self._algorithm_function = load_alignment_algorithm(self._algorithm)

        return self._algorithm_function(seq_1.sequence, seq_2.sequence,
                                        self._substitution_matrix, **kwargs)

    def print_aligned_sequences(self):

        # find longest name of target sequence and aligned sequences (for display purposes)
        max_name = max(len(self.target.name), max(len(x) for x in self._collection.names))

        f = '{:>%d}: {}' % max_name
        f_score = '{:>%d}: {} (Score: {})' % max_name

        # store the final string in a list so that everything is printed at the end in one go
        final_string = list()

        final_string.append(f.format(self.target.name, self.target.sequence))
        final_string.append('-'*len(f))

        for seq in self._collection.names:
            final_string.append(f_score.format(seq, self._aligned_collection[seq], self._alignment_scores[seq]))

        print(*final_string, sep='\n')

    @property
    def target_sequence(self):
        return self.target.sequence

    @property
    def aligned_sequences(self):
        return self._aligned_collection

    @property
    def score(self):
        return self._alignment_scores

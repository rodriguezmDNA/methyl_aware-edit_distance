from typing import Tuple
import numpy as np 

class MethylAwareDistance:

    def __init__(self, expected: str, observed:str,
                 read_1: bool = True,
                 match: int = 0, mismatch: int = 1, indel: int =1):

        self.expected = expected
        self.observed = observed
        self.read_1 = read_1

        self.size_observedatch = match
        self.size_observedismatch = mismatch
        self.indel =  indel

        self.expected = self.expected.upper()
        self.observed = self.observed.upper()

        self.size_expected = len(self.expected) 
        self.size_observed = len(self.observed)
        self.equal_length = self.size_expected == self.size_observed

        self.equivalent_bases = self.base_equivalency()

    def base_equivalency(self) -> Tuple[str,str]: 

        """
        Return equivalent bases during methyl converexpectedion (C>T or G>A) depending on the read (read 1 or read 2).

        Methylated Cs are generally converted to T during conversion.
        If conversion efficiency is low, methylated Cs can still be converted to T.
        
        
        -----> Read 1
        ================ Fragment ==============
                                <----- Read 2                                


        ## Read 1
        5'---mCG---3'
            
        Bad protection:  5'---TG---3'
        Good protection: 5'---CG---3'


        ## Read 2
    
        5'---CpG---3'
        3'---GmC---5'  || Read 2: bottom strand
        

        3'---GT---5'  || C>T in the bottom strand


        5'---CA---3'  || G>A conversion when normalized to 5'>3'
        3'---GT---5'  || 
        """

        return ('C','T') if self.read_1 else ('G','A')

    def hamming_methyl_aware(self) -> int:
        """
        Calculate the methylation-aware Hamming distance between two sequences.
        
        Due to the nature of cytosine conversion, the order of the sequences matters.
        """
        
        if not self.equal_length:
            raise ValueError("Error: lengths of sequences should be the same.")

        else:
            distance = 0
            
            for base_exp, base_obs in zip(self.expected, self.observed):
                if (base_exp != base_obs) and not (base_exp, base_obs) == self.equivalent_bases:
                    distance += self.size_observedismatch
            
            return distance


    def methyl_aware_base_score(self, s_i: str, t_j: str) -> int: 
        '''
        Return the score for a match or indel/mismatch accounting for methyl conversion
        '''
        return self.size_observedatch if (s_i == t_j) or (s_i, t_j) == self.equivalent_bases else self.size_observedismatch

    
    def levenshtein_methyl_aware(self) -> int:


        distance_matrix = np.zeros((self.size_expected+1,
                       self.size_observed+1),
                       dtype=int)

        distance_matrix[0,1:] = self.indel
        distance_matrix[0,:] = np.cumsum(distance_matrix[0,:])
        
        distance_matrix[1:,0] = self.indel
        distance_matrix[:,0] = np.cumsum(distance_matrix[:,0])
        
        for i in range(1, self.size_expected + 1):
            for j in range(1, self.size_observed + 1):

                match_mismatch_score = self.methyl_aware_base_score(self.expected[i - 1],
                                                                    self.observed[j - 1]
                                                                    )
                

                # if (expected[i - 1] == observed[j - 1]) or (expected[i - 1],observed[j - 1]) == equivalent_bases:
                #     dp[i][j] = dp[i - 1][j - 1] + match
                # else:    
                indelmm = [
                        (distance_matrix[i - 1][j] + self.indel), # indel in the expected sequence
                        (distance_matrix[i][j - 1] + self.indel), # indel in the observed sequence
                        (distance_matrix[i - 1][j - 1] + match_mismatch_score) # mismatch in the 
                        ]
                
                distance_matrix[i][j]  = indelmm[np.argmin(indelmm)] 

        return int(distance_matrix[-1,-1]) 


    def fill_distance_matrix(self, 
                             distance_matrix: np.ndarray,
                             index_observed: int
                             ):
        """
        Helper function for methyl-aware and space-efficient Levenshtein distance calculation.
        Takes in a 2-column matrix and updates the distance between the observed and expected sequences at a j point.
        Output is the updated matrix.
        """
        
        for index_expected in range(1, self.size_expected+1):

            ## Calculate match or mismatch score based on current characters from expected and observed sequences
            match_mismatch_score = self.methyl_aware_base_score( self.expected[index_expected - 1],
                                                                             self.observed[index_observed - 1]
                                                                             )
            

            # Calculate possible scores from different operations:
            scores = [
                distance_matrix[index_expected, 0] + self.indel,    # score from insertion/deletion in expected
                distance_matrix[index_expected-1, 1] + self.indel,  # score from insertion/deletion in observed
                distance_matrix[index_expected-1, 0] + match_mismatch_score  # score from match/mismatch base in observed and expected
            ]

            # Select the score that minimizes the number of edits
            distance_matrix[index_expected, 1] = min(scores)

        return distance_matrix

    def space_efficient_levenshtein_methyl_aware(self) -> int:

        """
        Function to calculate methyl-aware Levenshtein distance using a space-efficient implementation.
        For long sequences, filling a dynamic programming matrix is not efficient.

        Due to the nature of cytosine conversion, the order of the sequences matters.
        """


        # Initialize dynamic programming matrix. Two columns for space efficiency
        
        distance_matrix = np.zeros(( self.size_expected + 1, 2 ))
        distance_matrix[:,0]  = np.array( [0] + [i * self.indel for i in range(1, self.size_expected+1) ])
        distance_matrix[:,0] = distance_matrix[:,0].reshape(self.size_expected+1)
        
        # Iterate over elements of the second sequence
        for index_observed in range(1, self.size_observed + 1):
            distance_matrix[0, 1] = index_observed * self.indel
            distance_matrix = self.fill_distance_matrix(distance_matrix, index_observed)
            distance_matrix[:, 0] = distance_matrix[:, 1]
        
        return int(distance_matrix[-1,-1])
"""
class segregation is responsible for generating stress scores for Individual/unique variants of the pedigree
based upon the target pattern(inheritance) and observed pedigree pattern at a particular position. 
Please following along the example to understand the role of each of the functions below:

target_pattern abaaab (paired form ab aa ab)

chr:1    Pos:1626898    Ref:A      Alt: G  GT:GQ:GQX:DP:DPF:AD    1/1:48:48:17:0:0,17  0/1:137:137:23:2:9,14   1/1:78:78:27:1:0,27 

Considering the GT's for all the samples from the pedigree, the genotypes for the above position would be GGAGGG(genotype_list)
Next we process unique allele's from the genotype list, which in above case are A,G

For A:

The function standardized_genotype_vector_with_reference_to_a_particular_allele would create genotype vectors considering 'A' and the genotype list(G G A G G G)
for every match between 'A' and allele's from the genotype list a genotype_vector 'b' would be added to the genotype_vector list and 'a' for every miss match.
So the genotype_vector list for A vs (G G A G G G) would be 'aabaaa'. The genotype are then grouped in pairs and then sorted canonically so aabaaa becomes aa ab aa.

Next the target_test subroutine accepts genotype_vector list and the target_pattern list, checks to see if length's of the two list are equal and passes the two lists 
to stress_between_genotype_vectors function which further passes individual pairs from genotype_vector list and target_pattern to stress_between_two_genotypes.

eg: stress_between_genotype_vectors would pass 'aa' from the genotype_vector list and 'ab' from the target_pattern list to the stress_between_two_genotypes function. 
    next it would pass 'ab' from the genotype_vector list and  'aa' from the target_pattern list to the stress_between_two_genotypes function, and then 'aa' and 'ab'.

From the inputs of genotype_vector list and  target_pattern list stress_between_two_genotypes would come up with a segregation score for the queried allele, 
which in this case is'A'. 

All the above steps would be repeated for 'G' too.
 
Please note: If the target_test finds that the lengths of the genotype_vector list and the target_pattern list are not equal, it will skip this allele/position and report it in the error log. 

"""
from __future__ import print_function
import re
class segregation():
     
    def __init__(self):
        self.seg=""
    
    def stress_between_two_genotypes(self,observed_genotype,expected_genotype):
        stress=""
        observed_genotype=''.join(observed_genotype)
        expected_genotype=''.join(expected_genotype)
        observed_genotype.lower()
        expected_genotype.lower()
        sequencing_error_coefficient = 1
        if observed_genotype == expected_genotype:
            stress = 0
        else:
            if expected_genotype == "aa":
                if observed_genotype == "ab" or observed_genotype == "ba":
                    stress = 0.5
                elif observed_genotype == "bb":
                    stress = 1
                elif observed_genotype == "nn":
                    stress = 0.2
                elif observed_genotype == "an" or observed_genotype == "na":
                    stress = 0.1
                elif observed_genotype == "bn" or observed_genotype == "nb":
                    stress = 0.6
                elif observed_genotype == "aa":
                    quit() #case covered above
                else:
                    quit()
            elif expected_genotype == "ab" or expected_genotype == "ba":
                if observed_genotype == "ab" or observed_genotype == "ba":
                    stress = 0
                elif observed_genotype == "bb":
                    stress = 0.5
                elif observed_genotype == "nn":
                    stress = 0.2
                elif observed_genotype == "an" or observed_genotype == "na":
                    stress = 0.2
                elif observed_genotype == "bn" or observed_genotype == "nb":
                    stress = 0.1
                elif observed_genotype == "aa":
                    stress = 0.5
                else:
                    quit()

            elif expected_genotype == "bb":
                if observed_genotype == "ab" or observed_genotype == "ba":
                    stress = 0.5
                elif observed_genotype == "bb":
                    quit() #case covered above
                elif observed_genotype == "nn":
                    stress = 0.2
                elif observed_genotype == "an" or observed_genotype == "na":
                    stress = 0.6
                elif observed_genotype == "bn" or observed_genotype == "nb":
                    stress = 0.1
                elif observed_genotype == "aa":
                    stress = 1
                else:
                    quit()

            elif expected_genotype == "an" or expected_genotype == "na":
                if observed_genotype == "ab" or observed_genotype == "ba":
                    stress = 0
                elif observed_genotype == "bb":
                    stress = 0.5
                elif observed_genotype == "nn":
                    stress = 0.1
                elif observed_genotype == "an" or observed_genotype == "na":
                    stress = 0
                elif observed_genotype == "bn" or observed_genotype == "nb":
                    stress = 0.1
                elif observed_genotype == "aa":
                    stress = 0
                else:
                    quit()

            elif expected_genotype == "bn" or expected_genotype == "nb":
                if observed_genotype == "ab" or observed_genotype == "ba":
                    stress = 0
                elif observed_genotype == "bb":
                    stress = 0
                elif observed_genotype == "nn":
                    stress = 0.1
                elif observed_genotype == "an" or observed_genotype == "na":
                    stress = 0.1
                elif observed_genotype == "bn" or observed_genotype == "nb":
                    stress = 0
                elif observed_genotype == "aa":
                    stress = 0.5
                else:
                    quit()
            elif expected_genotype == "nn":
                stress = 0
            else:
                quit()
        return stress * sequencing_error_coefficient
    
 
    def stress_between_genotype_vectors(self,observed_vectors,target_pattern):
        stress=0
        observed_vectors=re.findall('..',observed_vectors)
        target_vectors=re.findall('..',target_pattern)
        myrange=len(observed_vectors)
        for x in range(0,myrange):
            stress += self.stress_between_two_genotypes(observed_vectors[x],target_vectors[x])
        return(stress)    

    def target_test(self,observed_vectors,target_pattern):
        if len(observed_vectors)!= len(target_pattern):
	    return 'NA'
        else:
            stress=self.stress_between_genotype_vectors(observed_vectors,target_pattern)
            return stress


    def standardized_genotype_vector_with_reference_to_a_particular_allele(self,allele_array,candidate_allele):
        canonical_vector=""
        allele_vector=""
        for vector in allele_array:
            if vector == candidate_allele:
                allele_vector +='b'
            else:
                allele_vector +='a'
        allele_vector=map(''.join, zip(*[iter(allele_vector)]*2))
        for allele in allele_vector:
            canonical_vector +=''.join(sorted(allele))
        return canonical_vector

#!/tools/bin/python
from __future__ import print_function 
from __future__ import division
from impala.dbapi import connect
from segregation import segregation
from collections import Counter
import datetime
import sys
import os.path
import ConfigParser
import getopt


#records timing
today = datetime.date.today()

#getopt read inputs from the commandline

targets=''
configfile=''
outdir=''

try:
    opts, args = getopt.getopt(sys.argv[1:],"t:o:c:",["targets=","outdir=","config=","help"])
except getopt.GetoptError:
    print ("One or more paramters are incorrect, please use pedigree_impala_analysis.py --help to know more")
    sys.exit(2)

for opt, arg in opts:
    if opt in ("--h","--help"):
        print ('Usage:pedigree_impala_analysis.py --targets <comma sep targets> --config <impala.config> --outdir <ab>')
        sys.exit()
    elif opt in ("--t", "--targets"):
        targets = arg
    elif opt in ("--c", "--config"):
        configfile=arg
    elif opt in ("--o","--outdir"):
        outdir=arg
    else:
        print ("Type python pedigree_impala_analysis.py --help to know more")
        sys.exit(2)



if not (targets and outdir and configfile):
    print("One or more required paramters are missing, please use pedigree_impala_analysis.py --help to know more")
    sys.exit(2)


##This allows reading and parsing of the configuration file

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


kaviar_cutoff=0
score_cutoff=0
stress_cutoff=0
##Reading of the configuration file to
##retrieve users input to the script

Config = ConfigParser.ConfigParser()
Config.read(configfile)

subject_id=ConfigSectionMap("Inputs")["inclusion_order"]
stress_cutoff=float(ConfigSectionMap("Inputs")["stress_cutoff"])
kaviar_cutoff=float(ConfigSectionMap("Inputs")["qaf_cutoff"])
score_cutoff=float(ConfigSectionMap("Inputs")["score_cutoff"])
restrict=ConfigSectionMap("Inputs")["restrict_to_chr"]
db_user=ConfigSectionMap("Inputs")["db_user"]

subject_list=subject_id.split(',')


##header_out add output headers to the output file
def header_out(pattern,writer):
    writer.write ('#Program_name: '+  __file__ + "\n")
    writer.write('#stress_cutoff: '+ str(stress_cutoff) + "\n")
    writer.write('#kaviar_cutoff: '+ str(kaviar_cutoff) + "\n")
    writer.write('#score_cutoff: ' + str(score_cutoff) + "\n")
    writer.write('#Inclusion_order: ' + str(subject_list) + "\n")
    writer.write('#target_pattern: ' + pattern + "\n")
    writer.write('#date:time: ' + str(datetime.datetime.now()) + "\n")
    writer.write("#Chromosome\tPosition\tReference\tCandidate\tGenotypes\tGenotype_vectors\tCandidate_stress\tAllele_count_in_pedigree\tMax_score\tMin_score\tAverage_score\tKaviar\tDANN\tCMS\tClinvar_sig\tGene_name\tOverall_candidate_score" + "\n")

#Creates dict/hash for all target patterns
target_inheritance=targets.split(",")
target_hash={}
for pattern in target_inheritance:
    completeName = os.path.join(outdir, 'candidate_pattern_'+pattern + '.txt')
    writer_pattern = open(completeName, 'w')
    target_hash[pattern]=writer_pattern
    header_out(pattern,writer_pattern)

##Database connection to Impala
try:
    conn=connect(host='glados19', port=21050,database=db_user)
    DB = conn.cursor()
    DB.execute("SELECT VERSION()")
    results=DB.fetchone()
    print ("Connection successful")
except:
    print("Connection failed, check if connection parameters are correct")
    sys.exit()

segregation_object=segregation()  #Object for the segregation class

"""
Vividict allows creation of perl like
hashes of hashes or multi key level hash
User would always have to point their hash/dict to Vividict
Eg: hash_dict=Vividict()

"""
class Vividict(dict):
    def __missing__(self, key):
        value= self[key] = type(self)()
        return value

"""
The compute_sequence_quality_score_statistics returns
max,min and average scores seen at a position.
"""
def compute_sequence_quality_score_statistics(score_array):
    max_score=max(score_array)
    min_score=min(score_array)
    num_list=[float(x) for x in score_array]
    avg_score=sum(num_list)/len(num_list)
    return (max_score,min_score,avg_score)

"""
QAF_score_modifier provides returns numerical
results based upon the allele frequenies that are
supplied. This score is essential in ranking the 
candidates or for calculating the candidate score
"""
def QAF_score_modifier(QAF): 
    if QAF == 0:
        return 1
    if QAF < 0.00001:
	return 1/0.99
    if QAF < 0.0001:
       return 1/0.98
    if QAF < 0.01:
        return 1/0.9
    if QAF < 0.05:
        return 1/0.75
    if QAF < 0.1:
        return 1/0.5
    if QAF < 0.15:
        return 1/0.3
    return 1

"""
This returns a numerical value based upon
what the score value is supplied.
The returned value is essential in calculating candidate score
"""
def quality_score_adjustment(max_score):
    score = max_score
    return_value=0
    if score >=50:
        return 0  #anything above 50 is good
    elif score > 35:
        return_value=(50-score)/15
        return return_value
    else:
        return 4

"""
This is similiar to above score modifier.
The motivation behind this function is to return low value
for highly pathogenic variants.
"""

def dann_score_modifier(dann_score):
    if dann_score >= 0.995:
       return 0.1
    elif dann_score >= 0.98 and dann_score < 0.995:
       return 0.5
    elif dann_score >= 0.93 and dann_score < 0.98:
       return 0.7
    else:
       return 1

"""
The overall candidate score is a rank given to a variant based 
upon inputs of genetic_stress, QAF,DANN, CMS etc.
Lower candidate score means that a variant is likely to be more pathogenic
"""
def overall_candidate_score(genetic_analysis_stress,queried_allele_frequency,dann_score,max_score,cms):
    QAF_score_modified=QAF_score_modifier(queried_allele_frequency)
    modified_dann=dann_score_modifier(dann_score)
    score = ((genetic_analysis_stress + 1) * (1 + QAF_score_modified) * (1 + cms) * modified_dann)
    adjusted_quality_score = quality_score_adjustment(max_score)
    overall_score = score * (1 + adjusted_quality_score)
    return(overall_score) 

"""
process_alleles is the heart of the pedigree_analysis program.
1)This re-creates family specific genotypes. 
2)Processes genotype quality
3)Creates dict/hash to store kaviar,dann,gene and cms annotation.
4)Annotates individual variants 
5)Calls the overall candidate score
6)Filters the output based upon user specified thresholds.
Currently,genotypes with GT ./. or 0/0 or missing from the impala database
are converted to homozygous reference calls.
"""
def process_alleles(list):
    kaviar={}
    quality={}
    dann={}
    clinvar={}
    subject=Vividict()
    GT={}
    chromosome=''
    position=''
    reference=''
    cms=''
    gene={}
    for line in list:
        chromosome=line[1]
	position=line[2]
        reference=line[3]
        subject[line[0]][line[4]]=1
        GT[line[0]]=line[5]
        try:
            if line[7] is not None:
                kaviar[line[4]]=line[7]
                if kaviar[line[4]] != line[7]:
	            kaviar[line[4]]+=line[7]
        except KeyError:
            kaviar[line[4]]=0
        
        if chromosome=='M':
            dann['A']=0.1
            dann['T']=0.1
            dann['G']=0.1
            dann['C']=0.1            
        elif chromosome=='M' and stress_cutoff >= 2.0:
            dann['A']=1
            dann['T']=1
            dann['G']=1
            dann['C']=1 
        else:
            dann['A']=line[8]
            dann['T']=line[9]
            dann['G']=line[10]
            dann['C']=line[11]
	
        if line[6] is None:
            quality[0]=1
	else:
	    quality[line[6]]=1
        clinvar[line[4]]=line[13] 
        cms=line[12]
        if line[14] is not None:
	    gene[line[14]]=1
    
    quality_list=[]
    
    for gqx in quality:
        quality_list.append(gqx)
    
    quality_array=compute_sequence_quality_score_statistics(quality_list)    
    max_score=quality_array[0]
    average_score="%.0f" % (quality_array[2])
    min_score=quality_array[1]    
      
    maximum=max(dann, key=dann.get)  
    dann['max']=dann[maximum]

    freq=0
    for value in kaviar:
        freq += kaviar[value]

    reference_freq=(1.0 - freq)    
    kaviar[reference]=reference_freq

    gene_string=' '.join(gene.keys())
    allele_string=''
    allele_array=[]
    for member in subject_list:
        if member in subject:
           if member in GT:
               if GT[member]=="1/1":
                   alt=subject[member].keys()
		   alt=''.join(alt)
                   allele_array.append(alt)
		   allele_array.append(alt)
		   allele_string += alt
		   allele_string += ' '
		   allele_string += alt
		   allele_string += ' '
               elif GT[member]=="0/1":
                   alt=subject[member].keys()
		   alt=''.join(alt)
		   allele_array.append(reference)
		   allele_array.append(alt)
		   allele_string += reference
		   allele_string += ' '
		   allele_string += alt
		   allele_string += ' '
	       elif GT[member]=="1/0":
                   alt=subject[member].keys()
                   alt=''.join(alt)
                   allele_array.append(alt)
                   allele_array.append(reference)
                   allele_string += alt
                   allele_string += ' '
                   allele_string += reference
                   allele_string += ' '
               elif GT[member]=="1/2" or GT[member]=="2/1":
                   alt = subject[member].keys()
		   for allele in alt:
		       allele_array.append(allele)
		       allele_string += allele
		       allele_string += ' '
		   
        else:
            allele_array.append(reference)
	    allele_array.append(reference)
	    allele_string += reference
	    allele_string += ' '
	    allele_string += reference
            allele_string += ' '
    allele_counter=Counter(allele_array)   
    Unique_array=set(allele_array)
    
    for candidate in Unique_array:
        for target_pattern in target_inheritance:
            genotype_vectors=segregation_object.standardized_genotype_vector_with_reference_to_a_particular_allele(allele_array,candidate)
            candidate_stress=segregation_object.target_test(genotype_vectors,target_pattern)
            if candidate_stress=='NA':
                continue
            kaviar_score=0
            dann_score=0
            cms_range=''
            clin_sig=''
            if candidate in kaviar:
                kaviar_score=kaviar[candidate]
            else:
                kaviar_score=0
     
            try:
                if dann[candidate] is not None:
                    dann_score="%4f" % dann[candidate]
            except KeyError: 
     	        dann_score=0
         
            if len(reference) > 1 or len(candidate) > 1:
                dann_score=dann['max']
                if dann_score:
                    dann_score="%4f" % (dann_score)
            elif candidate==reference:
                dann_score=dann['max']
                if dann_score:
                    dann_score="%4f" % (dann_score)
       
	    if candidate in clinvar:
	        if clinvar[candidate] is None:
                    clin_sig=0
                else:
	            clin_sig=clinvar[candidate]
            else:
	       clin_sig=0
        
            if cms is None:     
                cms_range=0
            else:
                cms_range=1
            
            output=[]   
	    overall_candidate_score_value=0                
            overall_candidate_score_value=overall_candidate_score(candidate_stress,kaviar_score,float(dann_score),max_score,cms_range)
            overall_candidate_score_value= float(overall_candidate_score_value)
            if kaviar_cutoff >= float(kaviar_score):
                if stress_cutoff >= candidate_stress:
                    if score_cutoff >=overall_candidate_score_value:
			if kaviar_score==0:
			    kaviar_score=0
 			elif kaviar_score < 0.0001:
		            kaviar_score="%2E" % (kaviar_score)
		            kaviar_score="%4f" %(float(kaviar_score))
			else:
			    kaviar_score="%4f" %(kaviar_score)
		        overall_candidate_score_value= "%4f" %(overall_candidate_score_value)
                        output=[chromosome,position,reference,candidate,allele_string,genotype_vectors,candidate_stress,allele_counter[candidate],max_score,min_score,average_score,kaviar_score,dann_score,cms_range,clin_sig,gene_string,overall_candidate_score_value]
	                for lines in output:
                            target_hash[target_pattern].write(str(lines) + "\t")
                        target_hash[target_pattern].write("\n")    


"""
This is reponsible for running chromosome sepcific left outer joins
between the temp tables and annotation tables.

"""
def query_impala(chr):
    mapper=Vividict()
    counter=0
    list=[]
    
    try:
        query=("select f.*, kav.allele_freq,dann.score_a as A,dann.score_t as T,"
           "dann.score_g as G,dann.score_c as C, cms.start,clin.clin_sig,ucsc.gene_name from "
          "(select fam.subject_id, fam.chrom, fam.pos,fam.ref,fam.alt,fam.gt,fam.gq"  
          " from %s.temp as fam "
          "where fam.chrom='%s') as f "
          "left outer join p7_ref_grch37.kaviar_isb as kav "
          "on f.chrom=kav.chrom "
          "and f.pos=kav.pos "
          "and f.ref=kav.ref "
   	  "and f.alt=kav.alt "
          "left outer join p7_ref_grch37.dann as dann "
          "on f.chrom=dann.chrom "
          "and f.pos=dann.pos "
          "left outer join p7_itmi.cms_gt1 as cms "
          "on concat('chr',f.chrom)=cms.chrom "
          "and f.pos >= cms.start and f.pos <= cms.stop "
          "left outer join p7_ref_grch37.clinvar as clin "
          "on f.chrom=clin.chrom "
          "and f.pos=clin.pos "
          "and f.ref=clin.ref "
          "and f.alt=clin.alt "
          "left outer join p7_ref_grch37.ucsc_genes as ucsc "
          "on f.chrom=ucsc.chrom "
          "and f.pos >= ucsc.txstart and f.pos <= ucsc.txend "          
          "order by f.pos" % (db_user,chr)) 
        DB.execute(query)
        print ("Finished running annotation joins on chrom " + str(chr),end="\n")

    except:
        print("Join failed, please check if impala is online and all tables are present",end="\n")
        print("Since joined has failed I am dropping the temp table",end="\n")
	sys.exit()

    for row in DB:
        subject_id = row[0]
        chrom = row[1]
        pos = row[2] 
        if mapper[chrom][pos]:
	    list.append(row)
            counter+=1
        else: 
            if len(list) > 0:
                process_alleles(list)
            list=[]
            mapper[chrom][pos]=1
            list.append(row)
    
    process_alleles(list) 


"""
Creates a temporary table so that annotation joins could be done
This is created in db_user database as temp
"""
def create_table():
    print ("Creating a temp table on database " + db_user,end="\n")
    
    query_string=''
    for string in subject_list:
        query_string += "'" + string + "'" + ","
    query_string=query_string[:-1]
    
    try:
        DB.execute(("DROP TABLE IF EXISTS %s.temp") %(db_user))
        query=("create table %s.temp as"
              "(select * from p7_platform.wgs_illumina_variant where subject_id IN (%s))" % (db_user,query_string))
        DB.execute(query) 
	print ("Temporary table created on Database " + db_user,end="\n" )
    except:
        print("Couldn't create a temporary table are you sure the database information is correct",end="\n")
        sys.exit()
     
    if restrict:
       restrict_analysis=restrict.split(',')
       for chr in restrict_analysis:
           query_impala(chr)
    else:
        chrArray=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M']
        for chr in chrArray:
            query_impala(chr)

    drop_table()

"""
Once the analysis done the temp table is dropped.
"""
def drop_table():
    try:
        query=("drop table %s.temp" %(db_user))
        DB.execute(query)
        print ("Analysis is complete deleting temporary table",end="\n")
    except:
        print ("Failed to delete temporary table",end="\n")

##Program starts

#create_table()


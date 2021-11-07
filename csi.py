# Co-Splicing Index (CSI) metric calculation

## Add how to use

## Add how is CSI calculated

## Separate input and execution script

# You should merge all the data you have into one csv sheet. But take note:

# 1. Column "expt_unit" refers to experimental unit, which is usually the gene name.
# To be specific, experiment unit refers to a gene/ cell line/ exons combination.
# CSI will be calculated with data within the experimental unit.
# However, sometimes a gene is treated with different SSOs give different 
# isoforms. In such cases, give DIFFERENT gene names to each experiment. 

# For instance, refer to BAG6 where one experiment showed 3 exons involved, and
# the other showed 7 exons involved. The number of exons alone will warrant each
# experiment to be classified as one unique experimental unit. Also, isoform A,
# B and C for each experiment refers to different exon combinations. Therefore,
# their data should not be mixed to calculate CSI. As long as their expt_unit
# name is different, you can ensure that their CSI calulation is based on their
# own experimental unit (e.g. BAG6_expt1, BAG6_expt2).

# 2. The same controls can be used for different SSO treatments as long as 
# they have the same gene/ cell line/ exons combination.
# (e.g. One expt_unit = 5 rows of control + 5 rows of SSO_A + 5 rows of SSO_B 
# is acceptable; the CSI will be calculated based on SSO_A against control and 
# SSO_B against control)

# ExptUnit Error: Number of isoforms or exons for each condition (cell line/ 
# treatment) in experimental unit is not the same

# 3. If you wish to calculate CSI without care for cell type, give them the same
# name under cell_line column.

# 4. Name E1, E2, E3... En for exons involved. If we have up to E10 and one gene 
# only has 4 exons involved, then type NA in the cells for the remaining 6 exons.
# Do not leave cells blank or use NA anywhere else

# Developer notes
# Check all ##check
##check subsetting of dataset after adding column for biological replicates is correct
##check no columns before expt_unit and no columns after relative_abundance in 
#   input or else all index splicing will fail
##check all index slicing
##dangers: 
    # [1] lstA = lstB
    #     change to lstA will change lstB and calculations will be wrong
    # [2] take care of types within dataframe so slicing works 

# replicates in datasheet should be named 1, 2, 3... consequtively starting from 1
# Fill all the cells, do not leave blank or '-' in cells, except that you can
# leave NA in the column for exons if the number of exons in this gene 
# is less than the number of exons in another gene (see example CSV)

# Edge cases
# edge1 :na, nan, Na, 
# edge2: same expt_unit has different, exon numbers
# edge3: RA does not sum to 1
# edge4: same expt_unit has same isoform name


# User-defined parameters
directory = None
file_name = None
control_name = None # Name of control under "treatment" column
plot_name = None
datasheet_name = None

# Modules
import sys
from user_inputs import *
import pandas as pd
from os import chdir
from os import mkdir
from os.path import join
from numpy import NaN
import plotly.express as px
from math import log
from collections import OrderedDict
from itertools import combinations
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from scipy.stats import levene
from scipy.stats import normaltest

# Given a list of [(1.0, a), (3.0, b), (2.0, c)...], where a, b, c are floats,
# return (a, b, c)
def sort_replicates(lst):
    d = {}
    final_lst = []
    
    # Generate dictionary
    for x, y in lst:
        d[x] = y
        
    # Make sorted dictionary into OrderedDict to preserve the order
    od = OrderedDict(sorted(d.items()))
    
    # Return list of values in sorted OrderedDict
    for key, value in od.items():
        final_lst += [value]
    
    return final_lst

def initiate(directory, file_name):
    
    ## Import data
    print("Importing data from specified directory")
    chdir(directory)
    df = pd.read_csv(file_name)

    ## Check if NA values in dataframe is in valid format
    print("Check dataframe for invalid NAs")
    if "Na" in df.values or "na" in df.values or "nA" in df.values:
        print("")
        print("Warning: 'Na', 'na' or 'nA' detected in dataframe.")
        print("Please replace all non-applicable exons to NA or nan.")
        print("")
        print("If this warning is irrelevant, ")
        print("(i.e. Na/ na/ nA is name of expt_unit/ cell_line/ treatment/ isoform/),")
        txt = input("press y and Enter to continue")
        if txt.upper().replace(" ","") == "Y" or txt.upper.replace(" ","") =="YES" or \
            txt.upper().replace(" ","") == "YE" or txt.upper().replace(" ","") == "YEAH":
            pass
        else:
            sys.exit("NA Nomenclature Error: Please rename datasheet cells.")
    else:
        print("NA check passed")
    
    ## Check if exons are in valid format
    print("Check if exons are in valid format")
    col = list(df.columns)[5:-2]
    for i in range(len(col)):
        # Check if E is present
        if "E" not in col[i]:
            print("")
            print("Error: Exon columns need to be named with consequtive E{num} e.g. E1 E2 E3")
            sys.exit("Exon Column Nomenclature Error: Please rename exon columns.")
        # Check if numbers are consequtive (i.e. E1, E2, E3...)
        if i+1 == int(col[i].replace("E","")):
            pass
        else:
            print("")
            print("Error: Exon columns are not in consequtive numbers, e.g. E1 E3...")
            sys.exit("Exon Column Error: Please ensure exon columns are consequtive.")
    print("Exon format check passed")
    
    ## Check for biological replicates
    print("Processing biological replicates")
    # Add all duplicates into another dataframe
    df_dup = df[df.duplicated(subset=list(df.columns[0:-2]), keep=False)]
    
    # Make a dataframe with empty relative_abundance and replicate_ID
    df_dup_empty = pd.DataFrame(df_dup.iloc[:,:-2].values, columns = df.iloc[:,:-2].columns).drop_duplicates(keep = 'first')
    
    df_dup_empty.insert(len(df_dup_empty.columns),'replicate_ID', NaN)
    df_dup_empty.insert(len(df_dup_empty.columns),'relative_abundance', NaN)
    
    # Add replicates column to both dataframes
    df_dup_empty.insert(len(df_dup_empty.columns),'replicates', NaN)
    df_dup.insert(len(df_dup.columns),'replicates', NaN)
    
    # Add replicates of relative_abundance values to df_dup_empty
    df_dup_empty.loc[:,'replicates'] = df_dup_empty.loc[:,'replicates'].astype('object')
    for index1, row1 in df_dup_empty.iterrows():
        for index2, row2 in df_dup.iterrows():
            # Remove NAs in list for list comparison
            # e.g. lst1 = [a, b, c, nan] 
            #      lst2 = [a, b, c, nan]
            # Then lst1 =/= lst2 because of the nan, so we remove nan 
            clean_lst1 = list(row1[:-3].dropna())
            clean_lst2 = list(row2[:-3].dropna())
            
            # If all values in row (except for "replicate_ID", "relative_abundance",
            # and "replicates") are the same between df_dup and df_dup_empty,
            # then add replicates from df_dup to df_dup_empty's "replicates"
            if clean_lst1 == clean_lst2:
                if type(df_dup_empty.loc[index1, "replicates"]) != list:
                    df_dup_empty.at[index1, "replicates"] = [ (row2["replicate_ID"] ,str(row2["relative_abundance"])) ]
                else:
                    df_dup_empty.at[index1, "replicates"] += [ (row2["replicate_ID"] , str(row2["relative_abundance"])) ]
            else:
                pass
                
    # Sort replicate list: replicates column now has a list of e.g.
    # [(4.0, d), (3.0, c), (2.0, b), (1.0, a)], which we need to sort 
    # and return[a, b, c, d] using function sort_replicates(lst)
    for index, row in df_dup_empty.iterrows():
        df_dup_empty.at[index, "replicates"] = sort_replicates(df_dup_empty.loc[index, "replicates"])
    
    # Calculate mean using the replicates and add the mean to relative_abundance
    for index1, row1 in df_dup_empty.iterrows():
        # Convert the strings in replicates to floats 
        df_dup_empty.at[index1, "replicates"] = [ float(x) for x in df_dup_empty.loc[index1, "replicates"] ]
        # Calculate mean using the replicates
        df_dup_empty.loc[index1, "relative_abundance"] = sum(df_dup_empty.loc[index1, "replicates"])/len(df_dup_empty.loc[index1, "replicates"])

    # Add duplicates as rows with relative_abundance = mean, and keep all the 
    # relative_abundance values as a tuple in df['replicates']
    df.loc[:,'replicates'] = NaN
        
    # Drop all duplicates in original df
    df = df.drop_duplicates(subset= list(df.columns[0:-3]), keep = False)
    
    # Add filled df_dup_empty to df
    df = df.append(df_dup_empty)
    print("Biological replicates processed")
    
    # Drop rows if nan in all of the cells of "expt_unit","cell_line","treatment","isoform" & "exons" 
    df = df.dropna(subset=["expt_unit","cell_line","treatment","isoform","exons"], how='all')
    return df

# Generate treatment_list for all treatments to be compared against control
def make_trmt_list(df):
    print("Generating list of all possible treatment conditions")
    treatment_list = []
    for EU in list(set(df["expt_unit"])):
        for cell_type in list(set(df.query(f"expt_unit=='{EU}'")["cell_line"])):
            for treatment in list(set(df.query(f"expt_unit=='{EU}' & cell_line=='{cell_type}'")["treatment"])):
                if treatment != control_name:
                    treatment_list += [(EU, cell_type, treatment)]
    print("All possible treatment conditions generated")
    return treatment_list

# Given dfT and dfC, exon 1 and exon2, presence1/2 (must be integer 1 or 0)
# returns (diff_RA, contributing_isoforms, total_isoforms, RA_cont_spliced_isoform)
##check input to presence1/2 must be integer
def calc_diffRA(df_treatment, df_control, exon1, exon2, presence1, presence2, dfCont):
    ### Get sum of replicates for treatment and control
    # e.g. For treatment dataframe,
    #              exon1    exon2    exon3    replicates
    # isoform A    1        1        1        (0.2, 0.3, 0.25)
    # isoform B    1        1        0        (0.3, 0.4, 0.3)
    # where replicates refer to relative abundance of each isoform obtained 
    # per biological replicate: (replicate1, replicate2, replicate3).
    # To get relative abundance of all isoforms with exon1-exon2 of 1-1,
    # we sum the column in "replicates" to get (0.5, 0.7, 0.55).
    # For treatment we obtained (0.5, 0.7, 0.55), and we do the same for
    # control dataframe and, say, we get (0.1, 0.2, 0.17).
    # We can then do pairwise comparison between (0.5, 0.7, 0.55) and 
    # (0.1, 0.2, 0.17) to get p-value.
    
    ## Calculations for treatment dataframe 
    sub_dfT = df_treatment.query(f'{exon1}=={presence1} & {exon2}=={presence2}')
    treatment_list = sub_dfT["replicates"]
    if type(treatment_list) == pd.core.series.Series:
        treatment_list = list(treatment_list)
        
    # Possible input to the conditionals below:
        # If there are no replicates, one isoform, treatment_list is nan
        # If there are no replicates, multiple isoforms, treatment_list is [nan, nan ...]
        # If there are replicates, one isoform, treatment_list is [0.3, 0.34, 0.31...]
        # If there are replicates, multiple isoforms, treatment_list is [ [0.3, 0.34, 0.31...], [0.3, 0.34, 0.31...] ...]
    
    try:
        if type(treatment_list) == list:
            # if treatment list is nested list e.g. [ [0.3, 0.2], [0.3, 0.3]...]
            if type(treatment_list[0]) == list:
                # if first item is NaN
                if pd.isna(treatment_list[0][0]):
                    treatment_list = NaN
                    treatment_sum = sum(sub_dfT["relative_abundance"])
                    
                # if there are no NaN, len of treatment list > 1 means there are replicates,
                elif len(treatment_list[0]) > 1:
                    raise ReplicateDetected      
                    
                else:
                    print("Error: Unexpected input detected in calc_diffRA")
                    sys.exit()
            # if first item is NaN
            if pd.isna(treatment_list[0]):
                treatment_list = NaN
                treatment_sum = sum(sub_dfT["relative_abundance"])
                
            # if there are no NaN, len of treatment list > 1 means there are replicates,
            elif len(treatment_list) > 1:
                raise ReplicateDetected      
                
            else:
                print("Error: Unexpected input detected in calc_diffRA")
                sys.exit()
        
        # if not list, then it must be nan, if not then sys.exit to prevent calc error
        elif pd.isna(treatment_list):
            treatment_list = NaN
            treatment_sum = sum(sub_dfT["relative_abundance"]) 
            
        else:
            print("Error: Unexpected input detected in calc_diffRA")
            sys.exit()
            
        # Only two possible inputs to the conditionals below:
            # Nested list of len 1: [[0.2, 0.3...]]
            # Nested list of len >1: [ [0.2, 0.3], [0.15, 0.33] ]
    except ReplicateDetected:
        # If there is only one isoform, we get nested list of len 1
        
        if len(treatment_list) == 1:
            treatment_list = treatment_list[0]
            
        # If more than 1 isoform, we get a nested list of len > 1
        #      Replicate     1    2
        # E.g. isoform A [ [0.2, 0.3], 
        #      isoform B   [0.15, 0.33] ]
        # We can zip to put values from same replicate in a list
        # then sum them: after zip = [ [0.2, 0.15], [0.3, 0.33] ]
        #                after sum = [ 0.35, 0.63 ]
        elif len(treatment_list) > 1:
            zip_treatment_list = list(zip(*treatment_list))
            treatment_list = [sum(x) for x in zip_treatment_list]
            
        # treatment_sum is the mean of treatment_list
        treatment_sum = sum(treatment_list)/len(treatment_list)
        
        
        
        
    ## Calculations for control dataframe 
    sub_dfC = df_control.query(f'{exon1}=={presence1} & {exon2}=={presence2}')
    control_list = list(sub_dfC["replicates"])
    if type(control_list) == pd.core.series.Series:
        control_list = list(control_list)    
        
    # Possible input to the conditionals below:
        # If there are no replicates, one isoform, control_list is nan
        # If there are no replicates, multiple isoforms, control_list is [nan, nan ...]
        # If there are replicates, one isoform, control_list is [0.3, 0.34, 0.31...]
        # If there are replicates, multiple isoforms, control_list is [ [0.3, 0.34, 0.31...], [0.3, 0.34, 0.31...] ...]
    
    try:
        if type(control_list) == list:
            # if control list is nested list e.g. [ [0.3, 0.2], [0.3, 0.3]...]
            if type(control_list[0]) == list:
                # if first item is NaN
                if pd.isna(control_list[0][0]):
                    control_list = NaN
                    control_sum = sum(sub_dfC["relative_abundance"])
                    
                # if there are no NaN, len of treatment list > 1 means there are replicates,
                elif len(control_list[0]) > 1:
                    raise ReplicateDetected      
                    
                else:
                    print("Error: Unexpected input detected in calc_diffRA")
                    sys.exit()
            # if first item is NaN, means there are no replicates
            if pd.isna(control_list[0]):
                control_list = NaN
                control_sum = sum(sub_dfC["relative_abundance"])
                
            # if there are no NaN, len of control list > 1 means there are replicates,
            elif len(control_list) > 1:
                raise ReplicateDetected      
                
            else:
                print("Error: Unexpected input detected in calc_diffRA")
                sys.exit()
        
        # if not list, then it must be nan, if not then sys.exit to prevent calc error
        elif pd.isna(control_list):
            control_list = NaN
            control_sum = sum(sub_dfC["relative_abundance"]) 
            
        else:
            print("Error: Unexpected input detected in calc_diffRA")
            sys.exit()
            
        # Only two possible inputs to the conditionals below:
            # Nested list of len 1: [[0.2, 0.3...]]
            # Nested list of len >1: [ [0.2, 0.3], [0.15, 0.33] ]
    except ReplicateDetected:
        
        # If there is only one isoform, thus we get a list of len 1
        if len(control_list) == 1:
            control_list = control_list[0]
            
        # If more than 1 isoform, we get a nested list of len > 1
        #      Replicate     1    2
        # E.g. isoform A [ [0.2, 0.3], 
        #      isoform B   [0.15, 0.33] ]
        # We can zip to put values from same replicate in a list
        # then sum them: after zip = [ [0.2, 0.15], [0.3, 0.33] ]
        #                after sum = [ 0.35, 0.63 ]
        elif len(control_list) > 1:
            zip_control_list = list(zip(*control_list))
            control_list = [sum(x) for x in zip_control_list]
            
        # control_sum is the mean of control_list
        control_sum = sum(control_list)/len(control_list)
    
    diff_RA = treatment_sum - control_sum
    
    
    ## contributing_isoforms: number of isoforms used in sum of exon2:presence2 (control)
    # The lower the number of non-zero isoforms in control, in treatment, the more 
    # confident that the isoforms in control are contributing to the target co-spliced form
    sub_dfC_nonzero = dfCont.query('relative_abundance != 0') 
    contributing_isoforms = len(sub_dfC_nonzero)
    
    # Reset to save memory
    sub_dfC = None
    sub_dfC_nonzero = None
    
    ## cospliced_isoforms: the number of co-spliced isoforms in treatment
    dfC_cospliced = sub_dfT.query('relative_abundance != 0')    
    cospliced_isoforms = len(dfC_cospliced)
    
    # Reset to save memory
    sub_dfT = None
    dfC_cospliced = None
    
    ## (depreciated) RA_cont_spliced_isoform: equals to sum of exon2:presence2 (control)
    # If the relative abundance of the splice isoform in control is zero
    # then we are confident that any increase is due to SSO. 
    # Otherwise, the presence of co-spliced isoform in treatment and control
    # suggests that the quantity in the former might not have changed since
    # we are measuring relative abundance (not absolute quantity)
    RA_cont_spliced_isoform = NaN
    
    return (diff_RA, contributing_isoforms, cospliced_isoforms, RA_cont_spliced_isoform,
            treatment_list, control_list)

class ReplicateDetected(Exception):
    pass

# Function to calculate "diff_relative_abundance", "splicing_type"
##check using df_sub instead of df other than in input and df.query
def process_diffRA(df, EU, cell_type, trmt, df_final):
    # Subset df by EU and cell_type
    df_sub = df.query(f'expt_unit=="{EU}" & cell_line=="{cell_type}"').reset_index(drop=True)    
    
    ## Label exon column cells by their actual exon names
    
    # Get list of exons and fill up empty exon columns with nan 
    exon_list = list(set(df_sub['exons']))[0].replace(' ', '').split(',')
    exon_list = ['E' + str(x) for x in exon_list]
    col_len = len(df_sub.columns[5:-3])   
    if len(exon_list) < col_len:
        exon_list += [NaN]* (col_len - len(exon_list))
    
    # Rename exon columns
    new_col_names = ["expt_unit", "cell_line", "treatment", "isoform", "exons"] + \
        exon_list + ["replicate_ID", "relative_abundance", "replicates"]
    df_sub.columns = new_col_names
    
    ### Trim exon columns to remove columns without any splicing activity
    
    # Check which exon column has no splicing activity and drop them
    col = list(df_sub.columns)[5:-3]
    col = [x for x in col if pd.isnull(x) != True]
    
    del_list = ["replicate_ID"]
    for i in col:
        # Check if exon column always have the same value or
        # if exon column values have na
        if df_sub[i].isnull().values.any() or len(set(df_sub[i])) == 1: 
            del_list += [i]    
    df_sub = df_sub.drop(del_list, axis=1)
    
    # Drop columns with all NaN values
    df_sub = df_sub.dropna(axis=1, how='all')
    try:
        df_sub['replicates']
    except KeyError:
        df_sub.insert(len(df_sub.columns),'replicates', NaN)
        
    ### Get treatment df_sub, reset index
    dfT = df_sub.query(f'expt_unit=="{EU}" & cell_line=="{cell_type}" \
                   & treatment=="{trmt}"').reset_index(drop=True)
    
    ### Get control df_sub, reset index
    dfC = df_sub.query(f'expt_unit=="{EU}" & cell_line=="{cell_type}" \
                       & treatment=="{control_name}"').reset_index(drop=True)
    
    # Check if control is present.
    if len(dfC.index.values) == 0:
        print(f"Error: Control for expt_unit: '{EU}', cell_line: '{cell_type}',  treatment: '{trmt}' is not found.")
        print("Please include the control or ensure that the control name in input datasheet matches control_name user parameter.")
        ##check Perhaps can improve by allowing programme to continue but giving NA for this {EU}/ {cell_type}/ {trmt}
        sys.exit("Missing control")
    
    ### Make list of all possible exon pairs by combination (order does not
    #   matter)
    exon_pair_list = list(df_sub.columns)[5:-2]
    combi = list(combinations(exon_pair_list, r=2))

    ### If all row values are the same except for RA, then we have biological 
    #   replicates. Find mean and use it to calculate RA_diff (for n >= 2)
    #   Use datapoints to conduct pairwise test and give p-value (for n >= 3)
    
    # Data structure of exon1, exon2 presence/ absence and the type of co-
    # splicing they can undergo
    possibilities = ( ('1', '1', "co-exclusion"), 
                       ('0', '0', "co-inclusion"),
                       ('1', '0', "swap"), 
                       ('0', '1', "swap") )
    
    # For exon1 and exon2 in combinations of all exons
    for e1, e2 in combi:
        # For each e1-e2, there are four possible types of presence1-presence2
        # in control, which if present, associates with a type of co-splicing
        for p1C, p2C, splice_type in possibilities:
            # Reset all variables
            E_stat = None
            E_p = None
            N_stat1 = None
            N_p1 = None
            N_stat2 = None
            N_p2 = None
            u = None
            p = None
            out = None
            dfCont = None
            dfTreat = None
            
            # Check if unspliced isoform is present (nonzero RA) in control 
            print(f"Examining exon {e1} ({p1C}) and exon {e2} ({p2C}) in control")
            dfCont = dfC.query(f'{e1}=={p1C} & {e2}=={p2C}')
            if sum(dfCont['relative_abundance']) != 0:
                ## Perform RA diff calculation between the co-spliced form of 
                ## treatment and control
                p1T = swap_dict[p1C]
                p2T = swap_dict[p2C]
                # Check if co-spliced isoform is present (nonzero RA) in treatment; 
                # if present, run calc_diffRA
                dfTreat = dfT.query(f'{e1}=={p1T} & {e2}=={p2T}')
                if sum(dfTreat['relative_abundance']) != 0:
                    out = calc_diffRA(dfT, dfC, e1, e2, p1T, p2T, dfCont)

                    ## Pairwise tests
                    # Possible inputs to the following conditionals:
                        # If there are replicates: we have the sum of isoforms:
                        #   e.g. [0.1, 0.1, 0.2]
                        # If there are no replicates: nan
                    try:
                        if type(out[4]) == list:

                            # If list has 1 or two items (not enough for pair-
                            # wise test)
                            if len(out[4]) <= 2:
                                p = 1
                                
                            # If list has 3 or more replicates (do pairwise test)
                            elif len(out[4]) >= 3:
                                raise ReplicateDetected
        
                        # When there are no replicates, out[4] = nan
                        elif pd.isna(out[4]):
                            p = 1
                            
                        else:
                            print("Error: unexpected input in process_diffRA")
                            print(f"The input is {out[4]}")
                            sys.exit()
                    
                    except ReplicateDetected:   
                    # Perform Mann-Whitney U test if sample size is small
                    # Non-parametric, unpaired test
                        if len(out[4]) >= 3 and len(out[4]) < 30:
                            u, p = mannwhitneyu(out[4], out[5], use_continuity=True, 
                                             alternative='two-sided')
                        
                        # When sample size exceeds 30, by Central Limit Theorem,
                        # we might be able perform t-test; 
                        # Check equal variance and normality assumption first
                        elif len(out[4]) >= 30:
                            # Check if variance is equal or not by Levene's test
                            E_stat, E_p = levene(out[4], out[5], center='median', proportiontocut=0.05)
                            N_stat1, N_p1 = normaltest(out[4], axis=0, nan_policy='propagate')
                            N_stat2, N_p2 = normaltest(out[5], axis=0, nan_policy='propagate')
                            # If equal variance and normal distribution, perform t-test
                            if E_p > 0.05 and N_stat1 > 0.05 and N_stat2 > 0.05:
                                u, p = ttest_ind(out[4], out[5],equal_var=True, nan_policy='propagate', 
                                                 alternative='two-sided')
                            # If unequal variance or not normal, perform Mann-Whitney
                            else:
                                u, p = mannwhitneyu(out[4], out[5], use_continuity=True, 
                                             alternative='two-sided')
                    # Update df_final
                    df_final.loc[len(df_final)] = [EU,         # expt_unit
                                                   trmt,       # treatment
                                                   cell_type,  # cell_line
                                                   e1, e2,     # control exon 1 & 2
                                                   splice_type,# splicing type
                                                   out[1],     # contributing isoforms in control
                                                   out[2],     # cospliced isoforms in treatment
                                                   p,          # p-value
                                                   out[0],     # RA diff
                                                   conf_calculator(out[1],out[2],p)] # confidence score

def conf_calculator(contributing_isoforms, cospliced_isoforms, p): 
    a = 1/(contributing_isoforms)   
    b = 1/(cospliced_isoforms)   
    c = -log(p) 
    if c > 50:
        c = 50
    return ((a+b) * 25) + (c * 10)

def plot(df_final, plot_name, directory):
    fig = px.scatter(df_final, 
                     x="confidence_score", 
                     y="diff_relative_abundance", 
                     color="expt_unit", 
                     hover_data=['index'],
                     title="Co-Splicing Index",
                     labels={
                         "confidence_score": "Confidence score",
                         "diff_relative_abundance": "Difference in relative abundance of co-spliced isoform (treatment - control)",
                         "expt_unit": "Gene or experimental unit"},
                     template='seaborn')
    
    fig.update_traces(marker=dict(size=12,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))
    
    # Create directory for plots
    try:
        mkdir(join(directory + "/output_plots"))
    except FileExistsError:
        pass
    
    # Write file to plots directory
    fig.write_html(f"{directory}/output_plots/{plot_name}.html")

##check; check if calculations here are correct!!!!

# Execute
if __name__ == '__main__':
    # Read data and preprocess data
    chdir(directory)
    df = initiate(directory, file_name)
    treatment_list = make_trmt_list(df)  
    print("")
    
    # Create dictionary to change 0 to 1 and vice versa
    swap_dict = {'1':'0', '0':'1'}
    
    # Create dataframe for storing output
    df_final = pd.DataFrame(columns = ["expt_unit", 
                                   "treatment", 
                                   "cell_line",
                                   "control_exon1", "control_exon2", 
                                   "splicing_type", 
                                   "contributing_isoforms_count",
                                   "cospliced_isoforms",
                                   "p-value",
                                   "diff_relative_abundance",                                    
                                   "confidence_score"])
    
    # Start CSI calculation on dataframe
    print("Initiating Co-Spliced Index (CSI) calculations")
    for x, y, z in treatment_list:
        print(f'Calculating for experiment unit "{x}", cell line "{y}" and treatment "{z}"')
        process_diffRA(df, x, y, z, df_final)
        print("Completed")
    
    # Set index
    df_final.reset_index(level=0, inplace=True)
    
    print("")
    print("All CSI calculations completed.")
    
    # Plot
    print("Generating plot")
    plot(df_final, plot_name, directory)
    
    ## Export datasheet as CSV
    print("Exporting datasheet")
    # Create directory for CSV
    try:
        mkdir(join(directory + "/output_datasheets"))
    except FileExistsError:
        pass
    
    df_final.to_csv(path_or_buf = f"{directory}/output_datasheets/{datasheet_name}.csv", index=False)
    print("All reports generated.")
    print("Programme complete.")    
    
    
    
    
    
    
    
    

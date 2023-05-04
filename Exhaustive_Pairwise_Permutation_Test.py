"""
EXHAUSTIVE PAIRWISE PERMUTATION TEST
(version 1.0)
by Angelo Chan

This is a program for conducting a series of two-group permutation tests.

The file being analyzed can contain data from multiple experiments. For each
experiment, a series of round-robin pairwise analyses will be performed. Every
group will be compared against all other groups.

The tests being performed can either be a frequentist pairwise permutation test,
(Possibly known by other names) or a standard deviation estimation permutation
test.

A frequentist pairwise permutation test calculates the differences between all
possible group ID shuffling permutations and finds the percentage of them in
which the difference was greater than that of using the original group IDs.

A standard deviation estimation permutation test calculates the average and the
standard deviation of the group mean differences for all shuffling permutations
of the group IDs, and calculates a p-value based on this.

It must be noted that these are not true p-values, and should not be taken as
such. Permutation tests are usually employed when the the conditions required
for more conventional statistical methods are not met. They allow users to at
least have some kind of objective number to use, which is better than nothing.



The output file will be a TSV. Each row will contain the results of the
comparison(s) between two groups in an experiment. The columns will be:
    - Experiment ID
    - Group ID 1 (first group in comparison)
    - Symbol signifying whether group 1 has a higher value, or group 2
    - Group ID 2 (second group in comparison)
    - (All comparison results)
    - (All extra information which the user specified should be kept)



USAGE:
    
    python27 Exhaustive_Pairwise_Permutation_Test.py <input_file>
            <{input_format}> <col_no_exp_id> <col_no_group_id> <col_nos_data>
            [-o <output_file>] [-t <{test_type}>] [-d <directional>]
            [-h <{header}> <{keep}>] [-k <col_nos_keep>]



MANDATORY:
    
    input_path
        
        The filepath of the input file.
    
    input_format
        
        The file format of the input file. Acceptable options are:
            tsv - Tab-separated values
            csv - Comma-separated values
            ssv - Space-separated values
    
    col_no_exp_id
        
        The column number of the column which contains the experiment IDs.
        
        Note that the column system uses the 1-index system. (The first column
        is column 1)
    
    col_no_group_id
        
        The column number of the column which contains the group IDs.
        
        Note that the column system uses the 1-index system. (The first column
        is column 1)
    
    col_nos_data
        
        A comma separated list of the columns containing the data to be
        analyzed.
        
        Note that the column system uses the 1-index system. (The first column
        is column 1)



OPTIONAL:
    
    output_file
        
        The filepath of the output file.
        If no path is specified, results will be directly outputted to the
        console.
    
    test_type
        
        (DEFAULT: F)
        
        The type of test to perform. Acceptable options are:
            F - Frequentist
            S - Standard Deviation
    
    directional
        
        (DEFAULT: N)
        
        Whether or not the test is a directional one.
        If yes, the test tests the probability of the larger group being larger
        than the smaller one by the original difference or more.
        If no, the test tests the probability of the two groups being further
        apart by the orginal difference or more.

    header
        
        (DEFAULT: Y)
        
        Whether or not there are column headers in the data file.

    keep
        
        (DEFAULT: Y)
        
        Whether or not to keep the column headers in the output file, if they
        are present in the input file.
    
    col_nos_keep
        
        A comma separated list of the columns containing the data to be kept and
        outputted, wholesale. If none are specified, no additional data will be
        kept.
        
        Note that the column system uses the 1-index system. (The first column
        is column 1)



EXAMPLES EXPLANATION:
    
    1:
    The experiment ID is in column 2. The group IDs are in column 3. We want to
    compare the data in column 1. The comparison being performed is a
    Frequentist one.
    
    2:
    The experiment ID is in column 1. The group IDs are in column 2. We want to
    compare the data in columns 10, 11, and 12. We want to keep the annotations
    in columns 3, 5, and 9. The file contains column headers which are not to be
    kept. An output file is specified.
    
    3:
    The experiment ID is in column 2. The group IDs are in column 3. We want to
    compare the data in column 1. The comparison being performed is a Standard
    Deviation one.

EXAMPLES:
    
    python27 Exhaustive_Pairwise_Permutation_Test.py example_data.tsv tsv 2 3 1
            -t F
    
    python27 Exhaustive_Pairwise_Permutation_Test.py example_data.tsv tsv 1 2
            10,11,12 -o output_file.tsv -h Y N -k 3,5,9
    
    python27 Exhaustive_Pairwise_Permutation_Test.py example_data.tsv tsv 2 3 1
            -t SD

USAGE:
    
    python27 Exhaustive_Pairwise_Permutation_Test.py <input_file>
            <{input_format}> <col_no_exp_id> <col_no_group_id> <col_nos_data>
            [-o <output_file>] [-t <{test_type}>] [-d <directional>]
            [-h <{header}> <{keep}>] [-k <col_nos_keep>]
"""

NAME = "Exhaustive_Pairwise_Permutation_Test.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Defaults #####################################################################

DEFAULT__test = 1 #FREQUENTIST
DEFAULT__directional = False
DEFAULT__header = True
DEFAULT__keep = True


# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from File_Reader import *
from Table_File_Reader import *
from Subgrouped_Table_File_Reader import *

from Simple_Permutator import *



# Enums ########################################################################

class TEST:
    FREQ=1
    SDEV=2



# Strings ######################################################################

STR__test_type = "\nERROR: Invalid test type:\n\t{s}"



STR__metrics = """
                    Average score: {A}

            Experiments processed: {B}
                 Groups processed: {C}
                  Lines processed: {D}

                Columns processed: {E}
    
    Average groups per experiment: {F}
     Average lines per experiment: {G}
    
          Average lines per group: {H}
"""

STR__report_begin = "\nRunning Exhaustive_Pairwise_Permutation_Test..."

STR__report_complete = "\nExhaustive_Pairwise_Permutation_Test successfully "\
        "finished."




# Lists ########################################################################

LIST__frequentist = ["F", "f", "FREQUENTIST", "Frequentist", "frequentist",
        "FREQ", "Freq", "freq"]
LIST__standard_d = ["SD", "Sd", "sd", "S", "s", "STANDARDDEVIATION",
        "StandardDeviation", "standarddeviation", "STANDARD_DEVIATION",
        "Standard_Deviation", "standard_deviation", "STANDARD", "Standard",
        "standard", "SDEV", "SDev", "sdev", "S_DEV", "S_Dev", "s_dev"]



# Dictionaries #################################################################

DICT__test = {}
for i in LIST__frequentist: DICT__test[i] = TEST.FREQ
for i in LIST__standard_d: DICT__test[i] = TEST.SDEV



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Exhaustive_Pairwise_Permutation_Test(path_in, delim, col_exp, col_grp,
            col_data, output_file, test_type, directional, header, keep,
            col_keep):
    """
    For each experiment, perform pairwise tests between 
    @path_in
            (str - filepath)
            The filepath of the input file.
    @delim
            (str)
            The delimiter to be used for the left table file. File formats and
            their corresponding delimiters are as follows:
                TSV - "\t" (tab character)
                CSV - ","  (comma character)
                SSV - " "  (whitespace character)
    @col_exp
            (int)
            The column number for the column containing the experiment IDs.
            (Uses a 0-index system.)
    @col_grp
            (int)
            The column number for the column containing the group IDs.
            (Uses a 0-index system.)
    @col_data
            (list<int>)
            A list of the column numbers for the columns which contain the
            values to be analyzed.
            (Uses a 0-index system.)
    @output_file
            (str - filepath)
            The filepath of the file where the output will be written into.
    @test_type
            (int - ENUM)
            An integer denoting what kind of test will be performed on the data.
            The options are as follows:
                1 - Frequentist
                2 - Standard Deviation
    @directional
            (bool)
            Whether or not the tests should be directional or not.
    @header
            (bool)
            Whether or not there are headers in the input files.
    @keep
            (bool)
            Whether or not to keep the headers, if they are present.
    @col_keep
            (list<int>)
            A list of the column numbers for the columns which contain the
            values to be kept as is. Note that only the values from the first
            row of data from each experiment will be kept.
            (Uses a 0-index system.)
    
    Exhaustive_Pairwise_Permutation_Test(path_in,delim,col_exp, col_grp,
            col_data,output_file,test_type,header,keep,col_keep)
    """
    PRINT.printP(STR__report_begin)
        
    PRINT.printP(STR__report_complete)
    
    # Reporting
    Report_Metrics(metrics)
    
    # Wrap up
    return 0



def Report_Metrics(metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @metrics
            (list<int>)
            A list of summary metrics for the data, including:
    
    Report_Metrics(list<int>(6)) -> None
    """
    return None
    


# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Consolidate_GE_Files(raw_command_line_input):
    """
    Parse the command line input and call the
    Exhaustive_Pairwise_Permutation_Test function with appropriate arguments if
    the command line input is valid.
    """
    PRINT.printP(STR__parsing_args)
    # Remove the runtime environment variable and program name from the inputs
    inputs = Strip_Non_Inputs(raw_command_line_input, NAME)
    
    # No inputs
    if not inputs:
        PRINT.printE(STR__no_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Help option
    if inputs[0] in LIST__help:
        print(HELP_DOC)
        return 0
    
    # Initial validation
    if len(inputs) < 5:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Setup mandatory inputs
    path_in = inputs.pop(0)
    input_format = inputs.pop(0)
    raw_col_exp = inputs.pop(0)
    raw_col_grp = inputs.pop(0)
    raw_col_data = inputs.pop(0)
    
    # Validate mandatory inputs
    valid = Validate_Read_Path(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    delim = Validate_File_Format(input_format)
    if not delim:
        printE(STR__invalid_file_format.format(io = "input", s = input_format))
        return 1
    col_exp = Validate_Int_Positive(raw_col_exp)
    if col_exp == -1:
        printE(STR__invalid_column.format(s = raw_col_exp))
        return 1
    col_exp -= 1
    col_grp = Validate_Int_Positive(raw_col_grp)
    if col_exp == -1:
        printE(STR__invalid_column.format(s = raw_col_grp))
        return 1
    col_grp -= 1
    col_data = Validate_List_Of_Ints_Positive(raw_col_data, ",")
    if not col_data:
        printE(STR__invalid_columns.format(s = raw_col_data, "commas"))
        return 1
    
    # Set up rest of the parsing
    output_file = None
    test_type = DEFAULT__test
    directional = DEFAULT__directional
    header = DEFAULT__header
    keep = DEFAULT__keep
    col_keep = []
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        try: # Following arguments
            if arg in ["-o", "-t", "-d", "-k"]:
                arg2 = inputs.pop(0)
            elif arg in ["-h"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except: # Redundant in current version
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        # Flag-dependent response
        if arg == "-o":
            path_out = arg2
        elif arg == "-t":
            if arg2 in DICT__test:
                test_type = DICT__test[arg2]
            else:
                PRINT.printE(STR__test_type.format(s = arg2))
                PRINT.printE(STR__use_help)
                return 1
        elif arg == "-d":
            directional = Validate_Bool(arg2)
            if directional == None:
                PRINT.printE(STR__invalid_arg_for_flag.format("-d"))
                PRINT.printE(STR__use_help)
                return 1
        elif arg == "-h":
            header = Validate_Bool(arg2)
            keep = Validate_Bool(arg3) 
            if (header == None) or (keep == None):
                PRINT.printE(STR__invalid_arg_for_flag.format("-h"))
                PRINT.printE(STR__use_help)
                return 1
            header = valid
        else: #arg == "-k"
            col_keep = Validate_List_Of_Ints_Positive(arg2, ",")
            if not col_keep:
                PRINT.printE(STR__invalid_columns.format(s=arg2, d="commas"))
                PRINT.printE(STR__use_help)
                return 1
            annotation_only = valid
    
    # Validate output path
    if path_out:
        valid_out = Validate_Write_Path__FILE(path_out)
        if valid_out == 2: return 0
        if valid_out == 3:
            PRINT.printE(STR__IO_error_write_forbid)
            return 1
        if valid_out == 4:
            PRINT.printE(STR__IO_error_write_unable)
            return 1

    # Run program
    exit_code = Exhaustive_Pairwise_Permutation_Test(path_in,delim,col_exp,
            col_grp,col_data,output_file,test_type,directional,header,keep,
            col_keep)
    
    # Safe exit
    if exit_state == 0: return 0
    else:
        PRINT.printE(STR__use_help)
        return 1



def Validate_Write_Path__FILE(filepath):
    """
    Validates the filepath of the output file.
    Return 0 if the filepath is writtable.
    Return 1 if the user decides to overwrite an existing file.
    Return 2 if the user declines to overwrite an existing file.
    Return 3 if the file exists and the program is set to forbid overwriting.
    Return 4 if the program is unable to write to the filepath specified.
    
    Validate_Write_Path(str) -> int
    """
    try:
        f = open(filepath, "U")
        f.close()
    except: # File does not exist. 
        try:
            f = open(filepath, "w")
            f.close()
            return 0 # File does not exist and it is possible to write
        except:
            return 4 # File does not exist but it is not possible to write
    # File exists
    if WRITE_PREVENT: return 3
    if WRITE_CONFIRM:
        confirm = raw_input(STR__overwrite_confirm.format(f = filepath))
        if confirm not in LIST__yes: return 2
    # User is not prevented from overwritting and may have chosen to overwrite
    try:
        f = open(filepath, "w")
        f.close()
        if WRITE_CONFIRM: return 1 # User has chosen to overwrite existing file
        return 0 # Overwriting existing file is possible
    except:
        return 4 # Unable to write to specified filepath



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Pairwise_Permutation_Test(sys.argv)



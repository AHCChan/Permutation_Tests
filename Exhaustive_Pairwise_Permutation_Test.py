HELP_DOC = """
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
    - Group ID 2 (second group in comparison)
    - (All comparison results)
        - Each result comprises two values and thus two columns.
        - The first column gives the probability value.
        - The second column gives the group ID of the higher value group.
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

CRUDE_Z_TEST = False



# Defaults #####################################################################

DEFAULT__test = 1 #FREQUENTIST
DEFAULT__directional = False
DEFAULT__header = True
DEFAULT__keep = True


# Imported Modules #############################################################

if not CRUDE_Z_TEST:
    import scipy.stats
else:
    from Crude_Z_Test import *



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

STR__use_help = "\nUse the -h option for help:\n\tpython "\
"Exhaustive_Pairwise_Permutation_Test.py -h"

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

LIST__headers = ["EXPERIMENT", "GROUP_1", "GROUP_2"]



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
            col_data, path_out, test_type, directional, header, keep,
            col_keep):
    """
    For each experiment, perform pairwise tests between the experimental groups.
    
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
    @path_out
            (str - filepath) OR
            (None)
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
    
    Exhaustive_Pairwise_Permutation_Test(path_in, delim, col_exp, col_grp,
            col_data, output_file, test_type, header, keep, col_keep)
    """
    PRINT.printP(STR__report_begin)
    
    # Setup - Metrics
    count_total = 0
    count_tests = 0
    count_exp = 0
    count_grp = 0
    count_line = 0
    
    # Setup - File I/O
    f = Subgrouped_Table_Reader()
    f.Set_New_Path(path_in)
    f.Set_Delimiter(delim)
    f.Set_Group_ID_Column_No(col_exp)
    if header:
        f.Set_Header_Params([1])
    f.Open()
    if path_out: o = open(path_out, "w")
    else: o = None
    
    # Setup - Others
    cols = [col_exp] + [col_grp] + col_data
    
    # Header
    if header and keep:
        # Process string
        header_str = f.Get_Header_Text()
        headers = header_str.split(delim)
        headers[-1] = headers[-1][:-1]
        # Build
        sb = delim.join(LIST__headers)
        sb += Build_Header_String(headers, col_exp, col_grp, col_data,
                col_keep, delim)
        # Write
        Controlled_Output(sb, o)
    
    # Main loop
    while not f.EOF:
        f.Read()
        raw = f.Get()
        # Annotations
        annotations = Build_String(raw[0], col_keep, delim)
        # Core
        data = Get_Data(raw, cols)
        results = Pairwise_Analyses(data, test_type, directional)
        total_score, total_tests, result_strs = results
        # Output
        for values in result_strs:
            sb = delim.join(values)
            sb += delim + annotations
            Controlled_Output(sb, o)
        # Metrics
        count_total += total_score
        count_tests += total_tests
        count_exp += 1
        count_grp += len(results)
        count_line += len(data)
    
    # Finish
    if path_out: o.close()
    f.Close()
    
    # Reporting
    Report_Metrics([count_total, count_tests, count_exp, count_grp, count_line,
            len(col_data)])
    
    # Wrap up
    return 0

def Get_Data(nested_lists, indexes):
    """
    Take a table of data (@nested_lists) and return only the relevant columns as
    specified by @indexes.
    
    @nested_lists
            (list<list<str>>)
            The data table, stored as a list of lists. Each nested list
            is one row of data.
    @indexes
            (list<int>)
            For each "row" of data, these are the indexes for the columns which
            contain the desired data.
            The first index is the column number for the experiment IDs, the
            second index is the column number for the group IDs, and the rest
            are the indexes for the data values to be analyzed.
    
    Get_Data(list<list<str>>, list<int>) -> list<[str, str, float/None...]>
    """
    result = []
    for sublist in nested_lists:
        temp = []
        temp.append(sublist[0])
        temp.append(sublist[1])
        for i in indexes[2:]:
            value = sublist[i]
            if value:
                value = float(value)
            else:
                value = None
            temp.append(value)
        result.append(temp)
    return result

def Pairwise_Analyses(data, test_type, directional):
    """
    Perform the relevant pairwise analyses and return the metrics of the
    resulting analysis and the strings to be output to the data file as a list
    of lists.
    
    @data
            (list<list>)
            A list of lists. Each sublist contains the data for one sample.
            The first value of the sublist is the experiment ID, a string.
            The second value is the group ID, also a string.
            The rest are data values as floats. If a value is missing, it will
            be a None.
    @test_type
            (int - ENUM)
            An integer denoting what kind of test will be performed on the data.
            The options are as follows:
                1 - Frequentist
                2 - Standard Deviation
    @directional
            (bool)
            Whether or not the tests should be directional or not.
    
    Pairwise_Analyses(<list<list>>, int, bool) -> [float, int, list<list<str>>]
    """
    # Setup - results
    results = []
    
    # Setup - reading
    length = len(data[0])
    length_ = length - 2
    range_ = range(length_)
    exp_ID = data[0][0]
    
    group_IDs = []
    for row in data:
        group_ID = row[1]
        group_IDs.append(group_ID)
    group_IDs = list(group_IDs)
    group_IDs = sorted(group_IDs)
    
    # Setup - metrics
    total_score = 0.0
    total_tests = 0
    
    # Subsets
    subsets = {}
    for group_ID in group_IDs:
        temp = []
        for i in range_:
            temp.append([]) # Necessary to ensure hardcopying
        subsets[group_ID] = temp
    for row in data:
        group_ID = row[1]
        for i in range(2, length):
            value = row[i]
            if value:
                subsets[group_ID][i-2].append(value)
    # Pairs
    pairs = Get_All_Pairs(group_IDs)
    
    for pair in pairs:
        # Unpack
        g1, g2 = pair
        
        # Setup
        row_result = [exp_ID, g1, g2]
        
        # All columns
        for i in range_:
            
            # From subset
            values_1 = subsets[g1][i]
            values_2 = subsets[g2][i]
            len_1 = len(values_1)
            len_2 = len(values_2)
            
            # Original
            avg_1 = sum(values_1)/len_1
            avg_2 = sum(values_2)/len_2
            if avg_1 > avg_2:
                larger = g1
                smaller = g2
                difference = avg_1 - avg_2
            else:
                larger = g2
                smaller = g1
                difference = avg_2 - avg_1
            
            # Setup permutations
            values = values_1 + values_2
            IDs = ([g1] * len_1) + ([g2] * len_2)
            permutations = Simple_Permutate(IDs)
            len_both = len_1 + len_2
            range_both = range(len_both)
            
            # Calculate all differences
            differences = []
            for permutation in permutations:
                total_g1 = 0
                total_g2 = 0
                for i in range_both:
                    value = values[i]
                    ID = permutation[i]
                    if ID == g1:
                        total_g1 += value
                    else:
                        total_g2 += value
                avg_1 = total_g1/len_1
                avg_2 = total_g2/len_2
                if larger == g1:
                    permutation_dif = avg_1 - avg_2
                else:
                    permutation_dif = avg_2 - avg_1
                differences.append(permutation_dif)

            # Calculate p-value
            p_value = Calculate_P_Value(difference, differences, test_type,
                    directional)
            total_score += p_value
            total_tests += 1
            
            # Stringbuilding
            row_result.append(str(p_value))
            row_result.append(larger)
            
        results.append(row_result)
    
    # Return
    return [total_score, total_tests, results]

def Calculate_P_Value(difference, differences, test_type, directional):
    """
    Calculate the p-value of a difference given a number of differences
    generated by permutation.
    
    @difference
            (float)
            The actual difference observed.
    @differences
            (list<float>)
            The differences generated by permutating the group IDs.
    @test_type
            (int - ENUM)
            An integer denoting what kind of test will be performed on the data.
            The options are as follows:
                1 - Frequentist
                2 - Standard Deviation
    @directional
            (bool)
            Whether or not the tests should be directional or not.
    
    Calculate_P_Value(float, list<float>, int) -> str
    """
    length = len(differences)
    if test_type == TEST.FREQ:
        count = 0
        for i in differences:
            if directional:
                if i >= difference:
                    count += 1
            else:
                if i >= abs(difference):
                    count += 1
        p = float(count)/length
        return p
    elif test_type == TEST.SDEV:
        mean = sum(differences)/length
        total = 0
        for i in differences:
            x = (i - mean) ** 2
            total += x
        sd = (total/length) ** 0.5
        z_score = difference/sd
        if directional:
            p = Flexible_Z_Test(z_score)
        else:
            p = Flexible_Z_Test(z_score) * 2
        return p
    return "NA"

def Flexible_Z_Test(z_score):
    """
    Convert a z-score into a p-value.
    
    Attempts to use the SciPy package, but can use an alternative if SciPy is
    not available.
    
    @z_score
            (float)
            Z-Score. The number of standard deviations a value is from the mean.
    
    Flexible_Z_Test(float) -> float
    """
    if CRUDE_Z_TEST:
        return SF(z_score)
    else:
        return scipy.stats.norm.sf(z_score)

def Build_Header_String(list_, col_exp, col_grp, col_data, col_keep, delim):
    """
    Build an output string for the header, according to the values and indexes
    given.
    
    @list_
            (list<str>)
            The raw values.
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
            The indexes of the columns containing data.
    @col_keep
            (list<int>)
            The indexes of the columns containing annotations.
    @delim
            (str)
            The delimiter used to separate the values.
    
    Build_String(list<str>, list<int>, str) -> str
    """
    headers_data = []
    for i in col_data:
        headers_data.append(list_[i])
    headers_keep = []
    for i in col_keep:
        headers_keep.append(list_[i])
    result = (list_[col_exp] + delim + list_[col_grp] + delim +
            (2*delim).join(headers_data) + delim + delim.join(headers_keep))
    return result

def Build_String(list_, indexes, delim):
    """
    Build an output string from the values in @list_, according to the indexes
    given. The values will be separated be @delim.
    
    @list_
            (list<str>)
            The raw values.
    @indexes
            (list<int>)
            The indexes of the values to be used.
    @delim
            (str)
            The delimiter used to separate the values.
    
    Build_String(list<str>, list<int>, str) -> str
    """
    headers = []
    for i in indexes:
        headers.append(list_[i])
    result = (delim).join(headers)
    return result     

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

def Controlled_Output(string, output_file):
    """
    Output the string to the output file, if an output file were specified.
    If no output file were specified, print the string to console.
    
    @string
            (str)
            The string to be outputted.
    @output_file
            (file) OR
            (None)
            The output file to which the string is to be written, or None.
    
    Controlled_Output(str, file) -> None
    Controlled_Output(str, None) -> None
    """
    if output_file:
        output_file.write(string + "\n")
    else:
        print(string)



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Pairwise_Permutation_Test(raw_command_line_input):
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
    delim = Validate_Table_File_Format(input_format)
    if not delim:
        PRINT.printE(STR__invalid_table_format.format(io = "input",
                s = input_format))
        return 1
    col_exp = Validate_Int_Positive(raw_col_exp)
    if col_exp == -1:
        PRINT.printE(STR__invalid_column.format(s = raw_col_exp))
        return 1
    col_exp -= 1
    col_grp = Validate_Int_Positive(raw_col_grp)
    if col_exp == -1:
        PRINT.printE(STR__invalid_column.format(s = raw_col_grp))
        return 1
    col_grp -= 1
    col_data = Validate_List_Of_Ints_Positive(raw_col_data, ",")
    if not col_data:
        PRINT.printE(STR__invalid_columns.format(s = raw_col_data,
                d = "commas"))
        return 1
    col_data = [i-1 for i in col_data]
    
    # Set up rest of the parsing
    path_out = None
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
        else: #arg == "-k"
            col_keep = Validate_List_Of_Ints_Positive(arg2, ",")
            col_keep = [i-1 for i in col_keep]
            if not col_keep:
                PRINT.printE(STR__invalid_columns.format(s=arg2, d="commas"))
                PRINT.printE(STR__use_help)
                return 1
    
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
    exit_state = Exhaustive_Pairwise_Permutation_Test(path_in,delim,col_exp,
            col_grp,col_data,path_out,test_type,directional,header,keep,
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



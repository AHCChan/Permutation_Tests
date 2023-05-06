
This ReadMe contains the following:

    Module description
    System requirements
    Dependencies
    Installation instructions



DESCRIPTION

A program for performing a series of pairwise permutation tests between all
groups in an experiment, for all experiments in a file.



REQUIREMENTS (SYSTEM)

This program runs in Python 2. It will not work properly with Python 3.
Please ensure you have Python 2 installed on your computer.
Please ensure you are using the correct version of Python to run this program.



REQUIREMENTS (DEPENDENCIES)

This program requires the following files from the following modules:

    File_Reader module: (https://github.com/AHCChan/File_Reader)
        File_Reader.py
        Table_File_Reader.py
        Subgrouped_Table_File_Reader.py

    Permutations module: (https://github.com/AHCChan/Permutator)
        Simple_Permutator.py

    Python_Command_Line_Tools module: (https://github.com/AHCChan/Python_Command_Line_Tools)
        _Command_Line_Parser.py
        _Controlled_Print.py



REQUIREMENTS (MULTIPLE CHOICE)

This program requires a module capable of performing Z-Tests. For the best
precision and accuracy, it is recommended that users install the SciPy and
NumPy modules, which can be installed in a number of ways. More information can
be found here:

    https://www.w3schools.com/python/numpy/numpy_intro.asp
    https://docs.scipy.org/doc/scipy/tutorial/general.html

Alternatively, for users who have trouble installing NumPy and SciPy, I have
created an alternative module:

    Crude_Alternatives module: (https://github.com/AHCChan/Crude_Alternatives)
        Crude_Z_Test.py

To use it, open your Exhaustive_Pairwise_Permutation_Test.py file in a text
editor or code editor.

Look for the line of code: 
    CRUDE_Z_TEST = False

Change it to:
    CRUDE_Z_TEST = True

Also ensure you actually download the module.



INSTRUCTIONS (SETUP)

You don't need to know how to use git to get this program, nor the libraries 
which this program needs to run.

You can simply open the files:
    File_Reader.py
    Table_File_Reader.py
    Subgrouped_Table_File_Reader.py
    Simple_Permutator.py
    Crude_Z_Test.py
    _Command_Line_Parser.py
    _Controlled_Print.py

... in your online browser and copy their contents onto a local file on your
computer with the same name. Ensure that these files are in the same folder.

Supporting files are available at:
    https://github.com/AHCChan/File_Reader
    https://github.com/AHCChan/Permutator
    https://github.com/AHCChan/Crude_Alternatives
    https://github.com/AHCChan/Python_Command_Line_Tools



INSTRUCTIONS (HOW TO USE)

You don't need to know how to use git to use this program. You can simply open
the Python files in your online browser and copy their contents into a local 
file on your computer with the same name.

To run the program, open up a command line window and enter one of the 
following commands (substituting the appropriate file paths) to bring up the 
relevant help doc:

    C:\Path\To\Python\python.exe C:\Path\To\The\File\Exhaustive_Pairwise_Permutation_Test.py -h

You should get a large wall of text explaining how to use this program, along
with some examples.

Alternatively, code in these files can be used as a module by other Python 
programs using standard import methods.



OTHER USEFUL TOOLS

The Table_To_Table.py tool (https://github.com/AHCChan/Table_To_Table) is a
useful general purpose tool for dealing with data in a table file format. (Such
as TSVs and CSVs) Consider using it in conjunction with the tools in this
package.



TESTING AND FEEDBACK

Feel free to contact Angelo Chan (angelo.hoi.chung.chan@gmail.com) if you have
any questions, feedback, or bugs to report.



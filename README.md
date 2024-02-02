# FastTree
Fast Tree Project CS4255 Algorithms for Sequence-based Bioinformatics

## Instructions
To run the algorithm make sure to have python installed on your pc. Download the main.py file and put it in a directory of your choice. Open a terminal or cmd here and navigate to the folder in which you have stored the main.py file. Now run the command:
```
python main.py "[file location]"
```
Replace [file location] with the location of your data file, eg c:\users\data\test-small.aln. Make sure to embed the file location in double quotation marks as shown in the example above. This is to prevent splitting up of the file location in different arguments in case the file location contains a space. if no data file is provided it will look in the same directory of the main.py file and search for a file named "fasttree-input.aln".

When you want to run the tests, make sure the tests.py file is located in the same directory as the main.py file. Then to run the test open a cmd and navigate to the folder both .py files are located. Finally, run the command:
```
python tests.py
```
This should return the test output.
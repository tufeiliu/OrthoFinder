# New feature note
diamond_cm_options_table = """
The following matrices are supported by DIAMOND, with the default being BLOSUM62.
Interpretation: int (gap open) / int (gap extend)

+----------+---------------------------------+-------------------------+
|  Matrix  |         Supported values        |  Default gap penalties  |
+----------+---------------------------------+-------------------------+
| BLOSUM45 | (10-13)/3; (12-16)/2; (16-19)/1 |           14/2          |
| BLOSUM50 | (9-13)/3; (12-16)/2; (15-19)/1  |           13/2          |
| BLOSUM62 | (6-11)/2; (9-13)/1              |           11/1          |
| BLOSUM80 | (6-9)/2; 13/2; 25/2; (9-11)/1   |           10/1          |
| BLOSUM90 | (6-9)/2; (9-11)/1               |           10/1          |
| PAM250   | (11-15)/3; (13-17)/2; (17-21)/1 |           14/2          |
| PAM70    | (6-8)/2; (9-11)/1               |           10/1          |
| PAM30    | (5-7)/2; (8-10)/1               |           9/1           |
+----------+---------------------------------+-------------------------+

For example, the default gap open and gap extend penalties for BLOSUM62 are 11 and 1, respectively
Apart from the default gap penalties, BLOSUM62 also support gap extend to be either 2 or 1,
with different gap open values allowed. 
For example, when gap extend is 2, gap open can be chosen between 6 and 11.

/// USAGE ///
Available scoring matrices used by DIAMOND:

1. `orthofinder -f ExampleData --matrix BLOSUM45`
When no gap penalties are specified, OrthoFinder will use the default gap penalties.

2. `orthofinder -f ExampleData --matrix BLOSUM62 --gapextend 2 --gapopen 10`
When specifying the gap penalties, gap extend penalty must be defined before the gap open penalty

2. `orthofinder -f ExampleData --matrix BLOSUM62 --gapextend 2`
If the gap open penalty is unavailable, OrthoFinder will use the largest gap open in allowed range defined by the provided gap extend penalty. 

Custom scoring matrices:

`orthofinder -f ExampleData -S diamond_custom --custom-matrix scoring-matrix-file.txt --gapopen 10 --gapextend 2`
When use a custom scoring matrix, the search program needs be changed to the custom version.
In the meantime, the gap penalties must be provided. 
There are no default gap penalties available in OrthoFinder for a custom scoring matrix.
The order of the gap extend and gap open penalties are unimportant in this case.
All entries in a custom scoring matrix must be integers.
"""
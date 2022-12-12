# Submission Attempts

- submission.csv - weighted average: 4,2,2,2,1,1,1,1,1 (0.603)
- submission2.csv - without weights: 1,1,1,1,1,1,1,1,1 (0.578)
- submission3.csv - based on model performance: 0.471,0.393,0.494,0.297,0.408,0.292,0.363,0.361,0.196 (0.58)
- submission4.csv - median (0.558)
- submission5.csv - geometric mean (0.55)
- submission6.csv - weighted average: 4,2,2,2,2,1,1,1,1 (0.598)
- submission7.csv - weighted average: 5,3,3,3,1,1,1,1,1 (0.605)
- submission8.csv - weighted average: 6,4,4,4,1,1,1,1,1 (0.605)
- submission9.csv - weighted average: 6,5,5,5,1,1,1,1,1 (0.604)
- submission10.csv - weighted average: 6,4,4,4,3,3,1,1,1 (0.6)
- submission11.csv - weighted average: 5,3,3,3,1,1,1,1 + remove blosum model (0.603)
- submission12.csv - weighted average: 5,3,3,3,1,1,1,1,1 + DEL 4,4,1,1 (0.605)
- submission13.csv - weighted average: 5,3,3,3,1,1,1,1,1 + DEL 10,10,1,1 (0.605)
- submission14.csv - weighted average: 5,3,3,3,1,1,1,1,1 + DEL 1,1,1,1 (0.604)
- submission15.csv - weighted average: 7,3,9,3,1,1,1,1,1 (0.598)
- submission16.csv - weighted average: 7,3,3,3,1,1,1,1,1 (0.603)
- submission17.csv - weighted average: 7,3,3,3,3,1,1,1,1 (0.599)
- submission18.csv - weighted average: 5,3,1,3,1,1,1,1,1 (0.6)
- submission19.csv - weighted average: 5,3,3,3,1,3,1,1,1 (0.605)
- submission20.csv - weighted average: 5,3,3,3,3,3,1,1,1 (0.596)
- submission21.csv - weighted average: 5,3,3,3,1,3,3,1,1 (0.599)
- submission22.csv - weighted average: 5,3,3,3,1,3,1,3,1 (0.6)
- submission23.csv - weighted average: 5,3,3,3,1,3,1,1,3 (0.603)
- submission24.csv - weighted average: 5,3,5,3,1,3,1,1,1 (0.603)
- submission25.csv - weighted average: 5,3,3,5,1,3,1,1,1 (0.604)
- submission26.csv - weighted average: 5,3,1,5,1,3,1,1,3 (0.598)
- submission27.csv - weighted average: 5,3,1,5,1,3,1,1,1 (0.6)
- submission28.csv - weighted average: 5,5,1,1,1,3,1,1,1 (0.587)


For weighted average, the 9 numbers correspond to the coefficients for the 9 models, in the following order:

1. rosetta
2. rmsd
3. thermonet
4. plddtdiff
5. sasaf
6. plddt
7. demask
8. ddg
9. blosum


For DEL, the 4 numbers correspond to the coefficients for the 4 models used in deletion mutations, in the following order:

1. plddt
2. plddtdiff
3. rmsd
4. sasaf
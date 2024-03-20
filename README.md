# osteometric-sorting

A GNU Octave function for osteometric sorting of human remains for bones of the upper and lower limbs.

The `sorting_BE` function utilizes diaphyseal cross-sectional geometric properties of the long bones for pair-matching. It can be used on femur, tibia, humerus, and ulna bones depending on the descriptives provided in the respective .mat files. The function can be called with one or three arguments as:

```
[sorted] = sorting_BE("descriptives.mat", "left_side.csv", "right_side.csv")

[sorted, stats] = sorting_BE("descriptives.mat", "left_side.csv", "right_side.csv")

[sorted, stats, unsorted] = sorting_BE("descriptives.mat", "left_side.csv", "right_side.csv")
```

The function prints three .csv files with the sorted elements (`sorted.csv`), the statistical information of the analysis, i.e. the number of total elements, pairs etc (`stats.csv`), and the plausible remaining pairs (`unsorted.csv`) according to the number of outputs.

The current version also supports descriptives in .csv files.


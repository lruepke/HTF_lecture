## Where is the cases and how to archive

All the cases are located in the same folder of a lecture, e.g. cases for Lecture 4 are in directory of `source/lectures/L04/cases/Jupp_Schultz`. The basic case will be compressed as a .zip file (called by `:download:` in the rst file) in the same directory, this is processed by `casePackage.sh` in the `cases` folder. The `casePackage.sh` is called in the `Makefile`.

**Note:** don't change file name of `casePackage.sh` for other cases, and keep the same folder name styles, e.g. `lectures/L0*`, `lectures/L0*/cases`.
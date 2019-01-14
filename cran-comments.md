## Test environments
* local OS install (MacOS Mojave version 10.14.2), R 3.5.2
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs.
* For OS, there were no NOTEs. 
* For win-builder, there were 3 NOTEs:

1. * checking CRAN incoming feasibility ... NOTE
Maintainer: 'Raymond Danner <dannerR@uncw.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  SongEvo (11:14, 12:26, 13:5)

2. ** running examples for arch 'i386' ... [40s] NOTE
Examples with CPU or elapsed time > 10s
          user system elapsed
par.sens 10.31   0.00   10.33
SongEvo   9.86   0.39   10.68
par.opt  10.05   0.05   10.11

3. ** running examples for arch 'x64' ... [56s] NOTE
Examples with CPU or elapsed time > 10s
          user system elapsed
SongEvo  15.94   0.67   16.77
par.opt  15.63   0.07   15.69
par.sens 11.61   0.01   11.62

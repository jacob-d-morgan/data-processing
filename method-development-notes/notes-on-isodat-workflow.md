# ISODAT Workflow

## Preparation
Before reprocessing files in ISODAT it is necessary to check that ISODAT will output the correct number of decimal places. Close all ISODAT programs and open Configurator. From the options tab, open the global settings window and set the number of digits for each of the value formats to the maximum number and click OK. This ensures that the intensities that the ISODAT reprocessing reports are precise enough for delta calculations.

## Reprocessing Files
This step is relatively simple. Once the '.did' and '.caf' result files are sorted chronologically and grouped by year they can be ISODAT Reprocessed. The catch is that ISODAT trips up if one of the result files contains no cycles. If there was an error in the peak center or pressure balance that caused the method to never perform any integrations then the file will not be reprocessed and the program will crash and likely close.

To avoid this, first sort the files by size and check that files with the smallest file size (typically ~68 kB) contain at least one cycle. It is also worth checking files with relatively unique sizes, in between the most common sizes, as they also often have no cycles. Place files with no cycles in a separate sub-directory of "Error-Causing Files" so that they are still present but do not intefere with the re-processing. Now, simply sort the files by creation date again so that they are reprocessed in the right order (youngest first, or resort and save the '.csv' file in Excel afterwards) and reprocess.

If ISODAT still trips up because you missed a file with no cycles, it is easy to identify the problematic file by opening the '.csv' file and seeing what the last file to be reprocessed file was. The next file is the problematic one. Remove it and retry.

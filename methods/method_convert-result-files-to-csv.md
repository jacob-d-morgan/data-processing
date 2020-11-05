# Method for Converting ISODAT Result Files to MATLAB-compatible Text Files 
This document describes the workflow for converting ISODAT result files to a '.csv' file that can be imported into MATLAB or other data analysis software.

## Make a Copy of the Results Folder
#### Use `robocopy` to make a duplicate of `C:\Thermo\Isodat NT\Global\User\Dual Inlet\Results` on a memory stick or an external drive. 
  
 - Open the command window and enter the command `robocopy "C:\Thermo\Isodat NT\Global\User\Dual Inlet\Results" "D:\Copy of Results" *.did *.caf /E /COPY:DAT /DCOPY:DAT /LOG:"D:\robocopy-log.txt" /TEE /NFL`. Substitute the path to your external drive for `D:\Copy of Results`. This copies ISODAT result files the results directory and all of its sub-directories to `D:\Copy of Results` while preserving the file creation date, which you will need later on. See the [robocopy documentation](https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/robocopy) for more details.
 - The folder is large so you will need plenty of storage space available on the drive.

## Rename the Files
#### Use PowerShell to rename the files, prepending the creation date and their relative path from the results folder.

 - Open Windows PowerShell ISE and navigate to the directory where the files are stored using the command line.
 - Open and run the PowerShell script.

## ISODAT Reprocess the Files
#### Use ISODAT to convert the files to '.csv' for importing to MATLAB.

 - Ensure that the configurator settings in ISODAT will report intensities to an appropriate number of decimal places.
 - Open ISODAT Workspace and navigate to the directory containing the PowerShell-processed files.
 - In each sub-directory:
     - Sort the files by size and inspect files with small or unusual file sizes to make sure they contain at least one cycle.
     - Place files that do not contain at least one cycle into a separate "Error-Causing Files" sub-directory.
     - Reorder the files by creation date
     - Select all the files, right click, and select re-process. Choose the appropriate template and set the output file format to '.csv'
     - Click OK. Make sure to check back frequently to ensure the program is still runnning.
     - If the program crashes, this is typically because it cannot reprocess a file containing less than one cycle. Open the '.csv' file and find the last file that was reprocessed. In workspace, place the following file in the "Error-Causing files" sub-directory and repeat the steps above.

# Workflow for Handling the ISODAT result files
This document describes the workflow I have developed for converting the ISODAT '.did' and '.caf' files, a.ka. "result files", to '.csv' files that can be passed to MATLAB.

## Approach Specifications
The workflow simply needs to convert the data in the proprietary ISODAT file formats to a text-based file format that can be read by MATLAB. There should be multiple '.csv' files (one per year) to keep the file sizes manageable and so that I can only import into MATLAB the data that I need for a given script.

Conversion of the result files to '.csv' can be achieved easily (albeit slowly), using the built-in ISODAT Reprocess functionality. I could therefore re-process the result files, in situ, on the mass spec PC, and then copy the resulting '.csv' files to my laptop. However, this is not ideal as the re-processing of files from even a single year takes many hours and could cause ISODAT to crash. This could mess up an analysis if a sequence were running, limiting me only to times when the MS is not being used. Because of these limitations, and in order to maintain a master copy of the results folder in case of emergency, it would be preferable to make a copy of the results folder that I can work with on a different machine.

In order to work with the result files on a different machine, I need to make a copy of `C:\Thermo\Isodat NT\Global\User\Dual Inlet\Results`. I then need to group all the result files by the year that they were created so that they can be ISODAT Reprocessed. This requires that (1) the copying of the results folder must preserve the file attributes (i.e. creation date), and (2) the grouping of the result files must preserve information about their original folder. I must preserve the creation date in order to group and order the files. The modification date is updated any time the result files are ISODAT Reprocessed, which would cause files to be misordered. I must preserve the folder that each file comes from as this is the only record of which ice core the file relates to.

## Approaches
### Preserving File Attributes
#### Unsuccessful Approaches
When copying a file to a new drive, such as a USB or external SSD, Windows sets the creation date to the current time. However, various methods are suggested online for copying files whilst preserving the file attributes. I tested many of them and found them to not work:
 - "Hold Ctrl": Holding the ctrl key whilst dragging and dropping file. Did not work.
 - "Send to Briefcase": Did not work.
 - "Zipping": Appears to work at first. Actually, the creation date of the unzipped file is set to the modification date.
 - "Batch File": Prepend creation date and relative path to the file name using a batch file, then copy as normal. Works, but is slow and crashes frequently so requires supervision.
 
#### Using Robocopy
Robocopy (robust file copy) is a built-in Microsoft command line tool, capable of copying large numbers of files in a way that preserves any number of their attributes, specified by flags. It is installed by default from Windows 7 onwards, but not on the version of XP that runs on the lab PCs. It was challenging to find a download link online as it seems to have been removed from the Windows Download Center. However, I was able to find an open source program ([Easy Robocopy](http://tribblesoft.com/easy-robocopy/)) that includes robocopy as part of its installation (N.B. The most recent version of Easy Robocopy does not include robocopy.exe. Instead download V1.0.9). This version of Robocopy (XP027) is portable, so can be installed on XP by copying the executable file to a new PC on a flash drive. I verified this using the lab laptop.

Robocopy allows me to make a full copy of `C:\Thermo\Isodat NT\Global\User\Dual Inlet\Results` on an external drive, preserving the creation dates of both the files and directories. This is done by executing the following line at the windows command prompt:

`robocopy "C:\Thermo\Isodat NT\Global\User\Dual Inlet\Results" "D:\Copy of Results" *.did *.caf /E /COPY:DAT /DCOPY:DAT /LOG:"D:\robocopy-log.txt" /TEE /NFL `
##### Syntax
 - `robocopy source destination files /options`
 - `"C:\Thermo\Isodat NT\Global\User\Dual Inlet\Results"`: Source Directory - copies files from this location
 - `"D:\Copy of Results"`: Destination Directory - copies files to this location
 - `*.did *.caf`: Files - copies files matching the filter(s). Wildcards capture all result files.
 - Options (see [documentation](https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/robocopy) for more details):
     - `/E` - Copies subdirectories, including empty subdirectories
     - `/COPY:DAT` - Copies file *D*ata, *A*ttributes, and *T*imestamps
     - `/DCOPY:DAT` - Copies directory *D*ata, *A*ttributes, and *T*imestamps
     - `/LOG` - Writes the command output to a log file
     - `/TEE` - Also writes the command output to the console
     - `/NFL` - Does not log each individual file copied, only the directories (No File List)
     
### Preserving Relative Path Information
Preserving information about the sub-directorys of `...Dual Inlet\Results` that a file is stored in is important as it is the only record of which ice core a file relates to. The information needs to be included in the '.csv' file that will be passed to MATLAB and so the easiest way to achieve this is to include the information in the file name, which is written into the '.csv' file during the ISODAT Reprocessing.

I can prepend the relative path of the results files using a relatively simple Windows PowerShell script, which also prepends the date created and groups the files into subdirectories by year and sorts them chronologically. The set of files is then ready to be handled in ISODAT.

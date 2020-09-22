    <#
    .SYNOPSIS
        Preprocesses ISODAT files from the Results folder, preparing them for ISODAT reprocessing.

    .DESCRIPTION
        Copies all files with a '.did' or '.caf' file extension from a directory and recursively
        from its subdirectories. Renames the copied files, prepending their creation date and
        their relative path from the starting directory and groups them into sub-directories based
        on the year the file was created. The present directory is used, unless a different
        location is specified using the -workingDir flag.

    .EXAMPLE
        PS C:\src\ps> preprocessResultFiles.ps1

        Description:
        Renames and groups the result files in the current directory.

    .EXAMPLE
        PS C:\src\ps> preprocessResultFiles.ps1 -workingDir C:\data\result-files

        Description:
        Renames and groups the result files in 'C:\data\result-files'.

    .NOTES
        Author: Jacob D. Morgan
        Contact: jacobmorgan@hotmail.co.uk
        Date: 17th September 2020
    #>

# Allow user to enter a path as an argument. Otherwise, use the current location
Param ([Parameter(Mandatory=$False)][string]$workingDir = (Get-Location))

# Find the Working Folder
$parentFolder = Split-Path $workingDir -Parent
$workingFolder = $workingDir.Substring($parentFolder.Length + 1)

# Find the files
$filesFound = Get-ChildItem -Path $workingDir -Recurse -Include *.did, *.caf

# Check there were some files found
If ($filesFound.Count -gt 0)
{

    # Print Number of Files Found
    Write-Host ("Found " + $filesFound.Count + " result files in " + $workingDir + " and its subfolders.")

    # Make a New Directory for the Files
    New-Item -Path ($workingDir + "\..\PowerShell-Processed-Files-from-" + $workingFolder) -ItemType Directory

    # Loop through them
    $idx = 1
    foreach ($resultFile in $filesFound)
    {
        # Get the Creation Date
        $creationDate = $resultFile.CreationTime | Get-Date -f "yyyy-MM-dd_HHmm"
        $creationDateUtc = $resultFile.CreationTimeUtc | Get-Date -f "yyyy-MM-dd_HHmm"
        $creationYear = $resultFile.CreationTime.Year.ToString()

        # Get the Relative Path from the Working Directory
        $absPath = $resultFile.DirectoryName
        if ($absPath.Length -gt $workingDir.Length + 1)
        {
            $relPath = $absPath.Substring($workingDir.Length + 1)
        }else {
            $relPath = ""
        }

        # Tidy up the Relative Path
        $relPathStr = $relPath.Replace(" ","-")
        $relPathStr = $relPathStr.Replace("\","_")
        $relPathStr = $relPathStr

        # Set the Target Directory
        $targetDir = ($workingDir + "\..\PowerShell-Processed-Files-from-" + $workingFolder +"\" + $creationYear)
        if (!(Test-Path $targetDir))
        {
            New-Item $targetDir -type Directory
        }

        # Set the New File Name
        $newFileName = $creationDate + "_" + $creationDateUtc + "utc_" + $relPathStr + "_" + $resultFile.Name

        # Copy and Rename the File
        Copy-Item $resultFile -Destination ($targetDir + "\" + $newFileName)
        Write-Host -Object ("Copying file '" + $resultFile.Name + "' from '" + $workingFolder + "\" + $relPath + "' (File " + $idx + " of " + $filesFound.Count + ")")

        $idx = $idx + 1;

    } # foreach resultFile
}else {
    Write-Host -Object ("No result files found in " + $workingDir + ".`n")
} # if filesFound
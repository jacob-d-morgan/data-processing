## First Attempt
I tried to install ISODAT 3.0 first. The installer would not run and gave me the error “This program requires a minimum resolution of 1024*768”. An internet search suggested changing the display scaling (display settings) from 150% to 125%. This fixed the problem.

The ISODAT installer can only install one of either Workspace (Office Installation) or Instrument Control and Acquisition (Instrument Installation). I had to install both in order to get Workspace to open. Prior to the instrument installation, Workspace would open and then immediately close itself.

After I got Workspace working, I installed ISODAT 2.0 as well. This caused problems because (I think) it installed in the same folder as ISODAT 3.0, presumably overwriting some of the installed files. I decided to reinstall them both in separate folders and see if that was better.

## Isodat 2.0
Sometimes the installer gets stuck right at the beginning. It is open in the taskbar but there is no window visible. I seemed to be able to fix this by opening task manager, showing the details of the process, and then ending the version handler task. This caused a window to appear and everything seemed to run smoothly after that.
Workspace opens fine after installing both the office and instrument programs and the service pack. Running the configurator tells workspace where to look for methods/results in the file browser. Seemingly no need to create any gas configs or anything else.

Unfortunately, ISODAT 2.0 doesn’t seem to be able to open the .caf files from 2013 and before. I get a sharing permission error and a second error message that suggests that the file is “not a valid ISODAT object”. A brief google search offers some suggestions that don’t work and ultimately suggests that the files are corrupted in some way that prevent them from being opened.  They also cannot be reprocessed, due to similar error.

## Isodat 3.0
The install has similar quirks to the 2.0 install. It is necessary to end the version control task in Task Manager after it trips up the installer. It is necessary to install instrument control/acquisition and workspace separately. In order to get workspace to function, it is necessary to open configurator once. It doesn’t seem to matter which options are selected, just that some options are set.

ISODAT 3.0 is able to open the files that were seemingly corrupted in 2.0. The version of 3.0 on the lab laptop opens the files without the permissions error, but does not display properly. There are weird tabs in strange places in the window, and portions of it are transparent, back to the desktop screen behind. The version of 3.0 I am able to install on my laptop (without the service pack) opens and displays the files normally as far as I can tell.

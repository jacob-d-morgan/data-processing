Plan of Attack:
---------------

- [x] Locate & Copy Files: Locate all the ISODAT result files and copy to my laptop
    - Done: Copied XP Results folder as of 2019-03-20

- [ ] Catalog Files: Catalog files by core

- [ ] Clean-Up Files: Identify and weed out unecessary results from misc tests or failed runs
    - Use file size?
    - Do in code instead?
    
- [ ] Design export template: Make an export template that is as versatile as possible.
    - Think about how I am going to sort samples and their blocks into the correct order
    
- [ ] Reprocess files: Reprocess all files from each core into one large excel file using ISODAT on lab laptop
        
- [ ] Import Raw Data Code: Write code that imports the raw delta values and as much sample metadata as possible from the reprocessed excel file
    - Include checks that identify potentially bad samples, maybe by checking the voltages, number of rows per block and blocks per sample?
    - I need to be able to identify the sample as one of: can vs can, PIS, CS, LJA, Ice, Unknown.
    
- [ ] Find Corrections Code: Write functions to extract PIS, CS (inc. tube), and LJA samples and processes them into a correction ID for each sample
    - Can make date-dependent correction here
    - Have to plot and check these very carefully, especially CS and LJA as they affect hundreds of ice samples
    
- [ ] Apply Corrections Code: Write functions to apply the PIS, CS (inc. tube), and LJA corrections

- [ ] Export Data Code: Write function to consolidate the raw data into a final dataset, sorted by bottom depth

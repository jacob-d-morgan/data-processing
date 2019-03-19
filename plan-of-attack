Plan of Attack:
---------------

[ ]  Locate all the ISODAT result files and copy to lab laptop

[ ]  Identify and weed out unecessary results from misc tests or failed runs
    - Use file size?
    - Do in code instead?
    
[ ]  Design an export template
    - Think about how I am going to sort samples and their blocks into the correct order
    
[ ]  Reprocess all files into one large excel file using ISODAT on lab laptop
        
[ ]  Write code that imports the raw delta values and as much sample metadata as possible from the reprocessed excel file
    - Include checks that identify potentially bad samples, maybe by checking the voltages, number of rows per block and blocks per sample?
    - I need to be able to identify the sample as one of: can vs can, PIS, CS, LJA, Ice, Unknown.
    
[ ]  Write functions to extract PIS, CS (inc. tube), and LJA samples and processes them into a correction ID for each sample
    - Can make date-dependent correction here
    - Have to plot and check these very carefully, especially CS and LJA as they affect hundreds of ice samples
    
[ ] Write functions to apply the PIS, CS (inc. tube), and LJA corrections

[ ] Write function to consolidate the raw data into a final dataset, sorted by bottom depth

%% Sample Type Timeline
% Plots a 'timeline' of the different categories of the imported data. This
% allows me to identify the different batches of CS, LJA, and which ice
% cores they are applied to.
% It also allows me to manually plot the date of significant changes to the
% MS, such as a refocussing, or significant adjustments to the method, such
% as the implementation of the O2-consumption correction.

sampleTypeCats = {'Std Can','Chem Slope','LJA','SPC','WDC','GISP','Misc. Test'};
sampleType = setcats(repmat(categorical(missing),length(cycle_metadata.filename),1),sampleTypeCats);

sampleInclusionCats = {'Included','Not Included','PIS'};
sampleInclusion = setcats(repmat(categorical(missing),length(cycle_metadata.filename),1),sampleInclusionCats);

sampleInclusion(ismember(cycle_metadata.msDatetime,aliquot_metadata.msDatetime))='Included';
sampleInclusion(~ismember(cycle_metadata.msDatetime,aliquot_metadata.msDatetime))='Not Included';
sampleInclusion(contains(cycle_metadata.ID1,'PIS'))='PIS';

isIncluded = sampleInclusion=='Included' | sampleInclusion=='PIS';

sampleType(contains(cycle_metadata.filename,'SPICE','IgnoreCase',true))='SPC'; % 'SPICE' is a folder name so this line will also capture all the LJA and Cans etc. in the folder
sampleType(contains(cycle_metadata.filename,'WDC','IgnoreCase',true))='WDC'; % 'WDC' is a folder name so this line will also capture all the LJA and Cans etc. in the folder
sampleType(contains(cycle_metadata.filename,'GISP','IgnoreCase',true))='GISP';

sampleType(contains(cycle_metadata.filename,'LJA','IgnoreCase',true))='LJA';
sampleType(contains(cycle_metadata.filename,'CS','IgnoreCase',true) | contains(cycle_metadata.filename,'Chem','IgnoreCase',true) )='Chem Slope';
sampleType(contains(cycle_metadata.filename,'Aair','IgnoreCase',true) | contains(cycle_metadata.filename,'Can','IgnoreCase',true))='Std Can';

sampleType(contains(cycle_metadata.filename,'Test','IgnoreCase',true))='Misc. Test';

%% Sample Type Timeline
% Plots a 'timeline' of the different categories of the imported data. This
% allows me to identify the different batches of CS, LJA, and which ice
% cores they are applied to.
% It also allows me to manually plot the date of significant changes to the
% MS, such as a refocussing, or significant adjustments to the method, such
% as the implementation of the O2-consumption correction.

%% Define and Assign Categories
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

%% Plot Timeline

figure; hold on;
plot(cycle_metadata.msDatetime(~isIncluded),zeros(sum(~isIncluded),1),'^k','MarkerFaceColor','k')

plot(cycle_metadata.msDatetime(isIncluded & isundefined(sampleType)),0.1*ones(sum(isIncluded & isundefined(sampleType)),1),'.','Color',lineCol(9),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & isundefined(sampleType)),0.1*ones(sum(~isIncluded & isundefined(sampleType)),1),'x','Color',lineCol(9)*0.7,'MarkerSize',5)

plot(cycle_metadata.msDatetime(isIncluded & sampleType=='Misc. Test'),0.2*ones(sum(isIncluded & sampleType=='Misc. Test'),1),'.','Color',lineCol(9)*0.7,'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='Misc. Test'),0.2*ones(sum(~isIncluded & sampleType=='Misc. Test'),1),'x','Color',lineCol(10),'MarkerSize',5)

plot(cycle_metadata.msDatetime(isIncluded & sampleType=='Std Can'),0.3*ones(sum(isIncluded & sampleType=='Std Can'),1),'.','Color',lineCol(3),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='Std Can'),0.3*ones(sum(~isIncluded & sampleType=='Std Can'),1),'x','Color',lineCol(3)*0.7,'MarkerSize',5)
plot(cycle_metadata.msDatetime(isIncluded & sampleType=='Chem Slope'),0.35*ones(sum(isIncluded & sampleType=='Chem Slope'),1),'.','Color',lineCol(8),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='Chem Slope'),0.35*ones(sum(~isIncluded & sampleType=='Chem Slope'),1),'x','Color',lineCol(8)*0.7,'MarkerSize',5)
plot(cycle_metadata.msDatetime(isIncluded & sampleType=='LJA'),0.4*ones(sum(isIncluded & sampleType=='LJA'),1),'.','Color',lineCol(6),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='LJA'),0.4*ones(sum(~isIncluded & sampleType=='LJA'),1),'x','Color',lineCol(6)*0.7,'MarkerSize',5)

plot(cycle_metadata.msDatetime(isIncluded & sampleType=='SPC'),0.5*ones(sum(isIncluded & sampleType=='SPC'),1),'.','Color',lineCol(2),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='SPC'),0.5*ones(sum(~isIncluded & sampleType=='SPC'),1),'x','Color',lineCol(2)*0.7,'MarkerSize',5)
plot(cycle_metadata.msDatetime(isIncluded & sampleType=='WDC'),0.55*ones(sum(isIncluded & sampleType=='WDC'),1),'.','Color',lineCol(1),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='WDC'),0.55*ones(sum(~isIncluded & sampleType=='WDC'),1),'x','Color',lineCol(1)*0.7,'MarkerSize',5)
plot(cycle_metadata.msDatetime(isIncluded & sampleType=='GISP'),0.6*ones(sum(isIncluded & sampleType=='GISP'),1),'.','Color',lineCol(7),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & sampleType=='GISP'),0.7*ones(sum(~isIncluded & sampleType=='GISP'),1),'x','Color',lineCol(7)*0.7,'MarkerSize',5)

ax=gca;
ax.YTick = [0 0.1 0.2 0.3 0.35 0.4 0.5 0.55 0.6];
ax.YTickLabel = {'Excluded Samples','Unkown','Misc. Test','Std. Can','Chem Slope','LJA','SPC','WDC','GISP'};


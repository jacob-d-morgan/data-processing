%% Test Aliquot Detection
% This script tests various methods of identifying aliquots from one
% another for reshaping.
% No one method is perfect, instead a combination of two or more of these
% approached has to be used.

% -------------------------------------------------------------------------

%% Using the Block Sequence Row
% A decrease in the value of the sequence row occurs only at the start of a
% new sequence, which almost guarantees a new sample.

iSeqRow = [-1; diff(block_metadata.sequenceRow(:,1))]<0;
idxSeqRow = find(iSeqRow);

%% Using the Block AS Inlet Value
% The first and last blocks of samples run from the SIO 8-port inlet have a
% value for this parameter, which tells the software which valves to
% operate when pumping away the previous sample and admitting the next.
% Can-vs-can (including PIS) blocks and blocks in the middle of a sample
% generally have a value of 'None'.

% Remove the 'None' blocks so that the first and last blocks of each sample
% are directly adjacent
iHasValue = block_metadata.ASInlet(:,1)~="None"; 
asInletHasValue = block_metadata.ASInlet(iHasValue,1);

% Create a comparison string of AS Inlet values, offset vertically by one
strToCompare = strings(size(asInletHasValue));
strToCompare(1) = "flag";
strToCompare(2:end) = asInletHasValue(1:end-1);

% Compare the strings to find the index where AS Inlet changes
iASInlet = false(size(block_metadata.msDatetime,1),1);
iASInlet(iHasValue) = ~strcmp(asInletHasValue,strToCompare);

% Missing values give a false positive, set them back to false
iASInlet(ismissing(block_metadata.ASInlet(:,1))) = false;

idxASInlet = find(iASInlet);

%% Using the Block ID1
% A change in ID1 indicates either a new sample OR a PIS block. This method
% identifies a lot more aliquots of only one block, the PIS blocks, than
% the other methods.

% Create a comparison string of ID1 values, offset vertically by one
strToCompare = strings(size(block_metadata.ID1(:,1)));
strToCompare(1,:) = "flag";
strToCompare(2:end,:) = block_metadata.ID1(1:end-1,1);

% Compare the strings to find the index where ID1 changes
iID1 = ~strcmp(block_metadata.ID1(:,1),strToCompare);

% Missing values give a false positive, set them back to false
iID1(ismissing(block_metadata.ID1(:,1))) = false;

idxID1 = find(iID1);

%% Using the Block Method
% The methods that pump out the inlet and introduce the next aliquot of
% sample gas are run only at the start of a sample. Define and look for
% these methods.

sampleStartMethods = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"]; % Add more methods here if necessary

iMethod = false(size(block_metadata.msDatetime,1),1); iMethod(1) = true;
for ii=1:length(sampleStartMethods)
    iMethod = iMethod | block_metadata.method(:,1)==sampleStartMethods(ii);
end

idxMethod = find(iMethod);

%% Using the Block Script
% The scripts that pump out the inlet and introduce the next aliquot of
% sample gas are run only at the start of a sample. Define and look for
% these scripts.
% N.B. for the 2016, 2017, and 2018 files, this gives the same results as
% using the block method. This may not be the case for other files if the
% name of the methods or scripts changed.

sampleStartScripts = [
    "dual inlet\acquisition64_cycle_1block_stds.isl"; 
    "dual inlet\acquisition64_cycle_1block_sa_std_delay.isl"; 
    "dual inlet\acquisition64_cycle_1block_sa_std.isl"]; % Add more scripts here if necessary

iScript = false(size(block_metadata.msDatetime,1),1); iScript(1) = true;
for ii=1:length(sampleStartScripts)
    iScript = iScript | block_metadata.scriptName(:,1)==sampleStartScripts(ii);
end

idxScript = find(iScript);


%% Plot the Different Methods
figure; hold on;
plot(iASInlet)
plot(iID1.*0.95)
plot(iMethod.*0.9)
plot(iSeqRow.*0.85)
legend('AS Inlet','ID1','Method','Seq. Row','Orientation','horizontal','Location','south')

figure
hold on;
plot(idxASInlet,1:length(idxASInlet))
plot(idxID1,1:length(idxID1))
plot(idxMethod,1:length(idxMethod))
plot(idxSeqRow,1:length(idxSeqRow))
legend('AS Inlet','ID1','Method','Seq. Row','Orientation','horizontal','Location','south')

%% Compare to Method & Sequence Row
% So far I have been using a combination of the method and sequence row
% approaches to detect aliquots. Compare the other two methods to this
% approach by finding the extra aliquots that they find.

extraASInlet = setdiff(idxASInlet,[idxMethod; idxSeqRow]);
extraID1 = setdiff(idxID1,[idxMethod; idxSeqRow]);
extraID1(block_metadata.ID1(extraID1,1)=='PIS') = [];

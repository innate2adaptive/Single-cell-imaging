function [DiffCells_OverTime, TotalCells_OverTime] = Sim_div_diff(trial,T,DD, delayall, delaydiff, isLogistic)
%this function simulates cell proliferation and differentiation in microwells over T hours.
%each simulated microwell starts with a different number of cells (1 - 10).
%the cells can divide or die at rates which were obtained by single cell
%experiments.
%the cells can differentiate either in a constant rate or in a rate that logistically depends on the instantaneous number of cells in a microwell.

%function inputs:
%trial - Number of simulates microwells. default is 100
%T     - maximal simulated hours. default is 96
%DD    - division destiny. a division density can be assigned to the cells which will limit the number of divisions. if it is 0 (default) then cells have unlimited divisions
%delayall   - delay the division time of all cells. default is 1 (no delay)
%delaydiff  - delay division time of differentiated cells. default is 1 (no delay)
%isLogistic - simulate either logistic differentiation (1, default) or 0 (constant differentiation)

%function outputs:
%TotalCells_OverTime - total number of cells in each time point for each microwell
%DiffCells_OverTime  - total number of differentiated cells in each time point for each microwell

%parameters uploaded from the 'parameter_mat':
%logistic differentiation parameters: maximum, Nc, rate
%tempbeta - distribution used to add noise to the simulation
%weil1d and weil2d - Weibull distributions used to draw first and subsequent division times.


load ('parameter_mat.mat');

%set default parameters
trial      = 100;
T          = 96;
DD         = 0;
delayall   = 1;
delaydiff  = 1;
isLogistic = 1;

divisiondest = makedist('Normal','mu',0,'sigma',1.8); %distribution used to add noise

TotalCells_OverTime= zeros(trial*10, T); %total cells over time
DiffCells_OverTime = zeros(trial*10, T); %differentiated cells over time
idx                = repelem(1:10,trial)'; %index for placing cells in the _time matrices

for q = 1:10 %number of initial cells
    
    if ~isLogistic
        diffrate = rateCon; %set constant diffrate. if fitting to a logistic function then diffrate is set at every time point
    end
    
    for w = 1:trial %Every well
        cells        = nan(10000,7); %Generate an empty array
        cells(1:q,1) = 1; %set this box to "cell" (for cell counts)
        cells(1:q,2) = 1; %are alive
        cells(1:q,3) = 0; %are differentiated
        cells(1:q,5) = round(random(weil1d,q,1)); %First division time drawn from a Weibull distribution based on the single cell data
        %col 1: cell Y/N
        %col 2: cell Alive/Dead
        %col 3: cell differentiated/un-differentiated
        %col 4: cell divisions remaining
        %col 5: cell division time
        
        if DD % optional to assign division destiny
            cells(1:q,4) = round(DD*random(tempbeta,q,1)); %assigning the cell's initial DD's
        else
            cells(1:q,4) = 100;  %without DD's - unlimited divisions
        end
        cells(1:q,6) = rand(q,1); %barcode each founder cell - provides the option to measure number of divisions of each clone
        
        for t = 0:T-1 %every hour
            celleff = random(divisiondest)+nansum(cells(:,2)); %get the effective cell number. this introduces noise to the simulation
            if celleff <= 0
                celleff = 0.000001; %make sure the effective cell number is never negative
            end
            if isLogistic
                diffrate=(maximum)*1/(1+exp(-rate*(celleff-Nc))); %set the differentiation rate for the well
            end
            
            liveref = cells(:,2)==1; %find every living cell
            DDleft  = ~isnan(cells(:,4)); %find cells with visions left
            cells(liveref&DDleft,5) = cells(liveref&DDleft,5)-1; %subtract one from division time
            
            liveDiv  = cells(:,5)==0; %find all the cells about to divide
            
            eventref = find(liveDiv&liveref&DDleft); %find all live cells with divisions left
            
            cells(cells(:,4)==0, 5) = nan; %if no more divisions left convert to nan
            cells(cells(:,4)==0, 4) = nan; %if no more divisions left convert to nan
            
            if(~isempty(eventref))
                cells(eventref(:),5) = round(random(weil2d,length(eventref(:)),1)*delayall); %create new division time drawn from a Weibull distribution. if delayall > 1 increase division time of the cell
                for n = 1:length(eventref(:,1)) %for every live cell that is about to divide
                    BC = cells(eventref(n),6); %get the cell's barcode
                    if sum(cells(:,6) == BC) == 1 %if the cell is in the first division -> there are no other cells with that barcode
                        deathrate = 0.15*(37/11); %death rate is multiplied by the relative rate per cell cycle (37 h for first division)
                    else
                        deathrate = 0.15; %death rate is multiplied by the relative rate per cell cycle (11 h for subsequent divisions)
                    end
                    if rand() <= deathrate %check if cell dies
                        cells(eventref(n),2) = 0;
                    else
                        if rand() <= diffrate %check if cell differentiates
                            cells(eventref(n),3) = 1; %if it does, set differentiation to true, redundancy doesn't matter
                            cells(eventref(n),5) = round(cells(eventref(n),5).*delaydiff); %if delaydiff > 1 increase division time of the differentiated cell
                        end
                        if cells(eventref(n),4) ~= 0 %if cell has divisions remaining
                            cells(eventref(n),4) = cells(eventref(n),4)-1; %subtract a division
                            cells(nansum(cells(:,1))+1,:)   = cells(eventref(n),:); %copy over the cell
                            cells(nansum(cells(:,1)),5)     = round(random(divisiondest)+cells(nansum(cells(:,1)),5));%add noise to the division time of sister cell
                            if cells(nansum(cells(:,1)),5) <= 0 %make sure the division time is never negative
                                cells(nansum(cells(:,1)),5) = 1;
                            end %of if
                        end %of if
                    end % of if
                end % on n
            end %of if
            stats(t+1,1) = nansum(cells(:,2)); %total living cells
            stats(t+1,2) = nansum(cells(find(cells(:,2)),3)); %total differentiated cells
        end % of t
        idx1 = find(idx==q); %get the locations in the matrix that are allocated for the 'q' cell number
        DiffCells_OverTime(idx1(w),:)   = stats(:,2)'; %total differentiated cells over time
        TotalCells_OverTime(idx1(w),:)  = stats(:,1)'; %total living cells over time
    end %of w
end % of q
end % of function
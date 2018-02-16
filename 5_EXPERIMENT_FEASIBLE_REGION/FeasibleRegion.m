clear; close all; clc;

runnable = true;   %#ok % Always set to true if runnable

% Load the MAC Variables
run('../Configuration_MAC');

%% Feasible Region Experiment
% thList = [3 4 5 6];        % List of the initial PSR (1./thList)
thList = [5 10/3];         % List of the initial PSR (1./thList)
% thList = [0.2 0.4 0.6 0.8];  % List of the initial PSR (1./thList)
% thList = [0.1 0.2 0.3 0.4];  % List of the initial PSR (1./thList)
conf.Niter = 10000;           % Keep it above 1000 iterations for meaningful results
conf.ColorList = {[155 155 155]./255 [0 0 0]./255};
% conf.ColorList = {'b' , 'g'};
FeasibleRegion_experiment(conf,thList);

%% Closing Clause
runnable = false;

%% Function Feasible Region
function FeasibleRegion_experiment(conf,thList)
    Niter = conf.Niter;
    alpha = conf.Alpha;
    ColorList = conf.ColorList;
    
    % Load PSR/PER Tables (PRX,SINR) -> PSR
    if(conf.ABSCfg == 0);      tablePSRName = strcat(conf.TablePath,'TABLE_PER_ABS0');
    elseif(conf.ABSCfg == 1)
        if(conf.MCS == 0);     tablePSRName = strcat(conf.TablePath,'TABLE_PER_ABS1');
        elseif(conf.MCS == 1); tablePSRName = strcat(conf.TablePath,'TABLE_MCS1_ABS1');
        elseif(conf.MCS == 2); tablePSRName = strcat(conf.TablePath,'TABLE_MCS2_ABS1');
        elseif(conf.MCS == 3); tablePSRName = strcat(conf.TablePath,'TABLE_MCS3_ABS1');
        end
    end
    load(tablePSRName,conf.varPSRList{:});
    conf.PRXLIST  = PRXLIST;
    conf.SINRLIST = SINRLIST;
    conf.PERLIST  = PERLIST./100;
    
	% Load MAC Tables (PSR) -> Average time to Transmit
    tableMACName = strcat(conf.TablePath,conf.TableName);
    load(tableMACName,conf.varMACList{:});
    conf.THPERLIST  = THPERLIST;
    conf.THPSDULIST = THPSDULIST;
    [~,idxAPtoGO]   = min(abs(repmat(alpha,size(ALPHAORG)) - ALPHAORG));
    [~,idxGOtoP2P]  = min(abs(repmat(1-alpha,size(ALPHAORG)) - ALPHAORG));
    conf.THTXTIMELIST_APtoGO  = THTXTIMELIST(idxAPtoGO,:);  %#ok % AP->GO
    conf.THTXTIMELIST_GOtoP2P = THTXTIMELIST(idxGOtoP2P,:);      % GO->P2P
    conf.THTXTIMELIST = THTXTIMELIST(end,:);                     % AP->P2P (alpha=1)
    
    % Feasible Region for PER and Throughput Oriented criterias
    x = (0.2:0.01:1);
    [X,Y] = meshgrid(x);
    Z2 = (1./X) + (1./Y);  % Throughput Oriented Criteria
%     Z3 = 1./(X.*Y);        % PSR Oriented criteria 

    str1 = 'PSR_{AP-WC} ='; str2 = '%';

    for idx = 1:length(thList)

        if(idx>1); margin = 2/100; % Percentage variation from the initial result
        else;      margin = 2/100; % Percentage variation from the initial result
        end

        conf.PSR_minth = (1/thList(idx)) * (1-margin);
        conf.PSR_maxth = (1/thList(idx)) * (1+margin);
        results = []; results_f = [];
        parfor iter = 1:Niter
            [~,~,~,~,~,~,outFeasExp] = EFi_function(conf);
            if ~isempty(outFeasExp.cand)
                results   = [results ; outFeasExp.cand];      
                results_f = [results_f ; outFeasExp.final];   
                for id = 1:size(outFeasExp.cand,1)
                    fprintf('%.4f\t%.4f\t%.4f\t%.2f\t%.2f\n',...
                       outFeasExp.cand(id,1), ...
                       outFeasExp.cand(id,2), ...
                       outFeasExp.cand(id,3), ...
                       outFeasExp.cand(id,4)*1e-3, ...
                       outFeasExp.cand(id,5)*1e-3);
                end
                for id = 1:size(outFeasExp.final,1)
                    fprintf('\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\n',...
                       outFeasExp.final(id,1), ...
                       outFeasExp.final(id,2), ...
                       outFeasExp.final(id,3), ...
                       outFeasExp.final(id,4)*1e-3, ...
                       outFeasExp.final(id,5)*1e-3);
                end 
            end
        end

        for f = 1:2
            figure(f); subplot(ceil(length(thList)/2),2,idx); hold on;
            str = strjoin({str1,num2str(100*(1/thList(idx))),str2});
            if(f==1)
                % Candidate-Final Comparrison
                sc(1) = scatter(results(:,3),results(:,2),25,'filled','o','MarkerEdgeColor',ColorList{1},'MarkerFaceColor',ColorList{1});
                sc(2) = scatter(results_f(:,3),results_f(:,2),25,'filled','s','MarkerEdgeColor',ColorList{2},'MarkerFaceColor',ColorList{2});
                title(str);
            elseif(f==2)
                % Throughput improvement for Candidates (surf)
                xlin = linspace(min(results(:,3)),max(results(:,2)),100);
                ylin = xlin;
                [T,W] = meshgrid(xlin,ylin);
                P1 = griddata(results(:,3),results(:,2),results(:,4),T,W,'cubic');
                P2 = griddata(results(:,3),results(:,2),results(:,5),T,W,'cubic');
                surf(T,W,P1,'FaceAlpha',0.5,'EdgeColor','none'); hold on;
                surf(T,W,P2,'FaceAlpha',0.5,'EdgeColor','none'); hold on;
                plot3(results(:,3),results(:,2),results(:,4),'.','MarkerSize',10,'Color','b');
                plot3(results(:,3),results(:,2),results(:,5),'.','MarkerSize',10,'Color','g');
                zlabel('Throughput (Mbps)');
                title(str);
                view([24.4 19.6]);
            end

            Z2_l = Z2 - thList(idx)*ones(size(X,1),size(X,2));
            C = contours(X,Y,Z2_l,[0 0]);
            xL = C(1,2:end); yL = C(2,2:end); zL = interp2(X,Y,Z2,xL,yL);
            zL = mean(results(:,4)).*ones(size(yL));
            h1 = line(xL,yL,zL,'Color','k','LineStyle','-.','LineWidth',2);
%             Z3_l = Z3 - thList(idx)*ones(size(X,1),size(X,2));
%             C = contours(X,Y,Z3_l,[0 0]);
%             xL = C(1, 2:end); yL = C(2, 2:end); zL = interp2(X,Y,Z3,xL,yL);
%             zL = mean(results(:,4)).*ones(size(yL));
%             h2 = line(xL,yL,zL,'Color','k','LineStyle',':','LineWidth',2);
            xlabel('PSR_{GO-WC}');
            ylabel('PSR_{AP-GO}');
            if(f==1 && idx == 2)
                legend([sc h1],{'Relay candidates','Relay final','Relay candidacy criteria (Eq. (3))'})
            end
            axis tight; ylim([(conf.PSR_th-0.15) 1]); grid minor;
        end
    end
end
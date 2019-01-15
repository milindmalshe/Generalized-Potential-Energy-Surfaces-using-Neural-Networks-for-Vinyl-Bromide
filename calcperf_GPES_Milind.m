

function [perf,Ex_Milind,Vhat_Milind]=calcperf_GPES_Milind(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
                                net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
                                net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr,...
                                V, Q, clust_size, type, ...
        r_C1C2, r_C1H3, r_C1H4, r_C1H5, r_C1Br6, r_C2H3, r_C2H4, r_C2H5, r_C2Br6, r_H3H4, r_H3H5, r_H3Br6, r_H4H5, r_H4Br6, r_H5Br6)


%%%%%$$$$$$$$$$@@@@@
perf=0.0;

fr_H1H2= 0; fr_H1Br= 0; fr_H1H2= 0; ftheta_H1H2Br= 0;

fc=1;%sim(net_fc,rij);

% % for iQ=1:1:Q
%     fr_H1H2= sim(net_fr_HH,r_H1H2);
%     fr_H1Br= sim(net_fr_HBr,r_H1Br);
%     fr_H2Br= sim(net_fr_HBr,r_H2Br);
%     
%     ftheta_H1H2Br= sim(net_ftheta_HHBr,[r_H1H2; r_H1Br; r_H2Br]);
%     
% % end
% Vhat_Milind= (fr_H1H2+fr_H1Br+fr_H2Br+ftheta_H1H2Br);


fr_all= sim(net_fr_CC,r_C1C2)+ sim(net_fr_HC, r_C1H3)+ sim(net_fr_HC,r_C1H4)+ sim(net_fr_HC,r_C1H5)+ sim(net_fr_CBr,r_C1Br6)+...
    sim(net_fr_HC,r_C2H3)+ sim(net_fr_HC,r_C2H4)+ sim(net_fr_HC,r_C2H5)+ sim(net_fr_CBr,r_C2Br6)+...
    sim(net_fr_HH,r_H3H4)+ sim(net_fr_HH,r_H3H5)+ sim(net_fr_HBr,r_H3Br6)+...
    sim(net_fr_HH,r_H4H5)+ sim(net_fr_HBr,r_H4Br6)+...
    sim(net_fr_HBr,r_H5Br6);

% ftheta_all= sim(net_ftheta_HCC,[r_C1C2; r_C1H3; r_C2H3])+ sim(net_ftheta_HCC,[r_C1C2; r_C1H4; r_C2H4])+ ...
%     sim(net_ftheta_HCC,[r_C1C2; r_C1H5; r_C2H5])+ sim(net_ftheta_HCBr,[r_C1C2; r_C1Br6; r_C2Br6])+ ...
%     sim(net_ftheta_HHC,[r_C1H3; r_C1H4; r_H3H4])+ sim(net_ftheta_HHC,[r_C1H3; r_C1H5; r_H3H5])+ sim(net_ftheta_HCBr,[r_C1H3; r_C1Br6; r_H3Br6])+ ...
%     sim(net_ftheta_HHC,[r_C1H4; r_C1H5; r_H4H5])+ sim(net_ftheta_HCBr,[r_C1H4; r_C1Br6; r_H4Br6])+ ...
%     sim(net_ftheta_HCBr,[r_C1H5; r_C1Br6; r_H5Br6])+ ...
%     sim(net_ftheta_HHC,[r_C2H3; r_C2H4; r_H3H4])+ sim(net_ftheta_HHC,[r_C2H3; r_C2H5; r_H3H5])+ sim(net_ftheta_HCBr,[r_C2H3; r_C2Br6; r_H3Br6])+...
%     sim(net_ftheta_HHC,[r_C2H4; r_C2H5; r_H4H5])+ sim(net_ftheta_HCBr,[r_C2H4; r_C2Br6; r_H4Br6])+...
%     sim(net_ftheta_HHH,[r_H3H4; r_H3H5; r_H4H5])+ sim(net_ftheta_HHBr,[r_H3H4; r_H3Br6; r_H4Br6])+...
%     sim(net_ftheta_HHBr,[r_H3H5; r_H3Br6; r_H5Br6])+...
%     sim(net_ftheta_HHBr,[r_H4H5; r_H4Br6; r_H5Br6]);


% fphi_all= sim(net_fphi_HHCC,[r_C1C2; r_C1H3 ; r_C1H4; r_C2H3; r_C2H4; r_H3H4])+ sim(net_fphi_HHCC,[r_C1C2; r_C1H3 ; r_C1H5; r_C2H3; r_C2H5; r_H3H5])+...
%     sim(net_fphi_HHCBr,[r_C1C2; r_C1H3 ; r_C1Br6; r_C2H3; r_C2Br6; r_H3Br6])+ ...
%     sim(net_fphi_HHCC,[r_C1C2; r_C1H4 ; r_C1H5; r_C2H4; r_C2H5; r_H4H5])+ sim(net_fphi_HCCBr,[r_C1C2; r_C1H4 ; r_C1Br6; r_C2H4; r_C2Br6; r_H4Br6])+...
%     sim(net_fphi_HCCBr,[r_C1C2; r_C1H5 ; r_C1Br6; r_C2H5; r_C2Br6; r_H5Br6])+...
%     sim(net_fphi_HHHC,[r_C1H3; r_C1H4 ; r_C1H5; r_H3H4; r_H3H5; r_H4H5])+ sim(net_fphi_HHCBr,[r_C1H3; r_C1H4 ; r_C1Br6; r_H3H4; r_H3Br6; r_H4Br6])+...
%     sim(net_fphi_HHCBr,[r_C1H3; r_C1H5 ; r_C1Br6; r_H3H5; r_H3Br6; r_H5Br6])+...
%     sim(net_fphi_HHCBr,[r_C1H4; r_C1H5 ; r_C1Br6; r_H4H5; r_H4Br6; r_H5Br6])+...
%     sim(net_fphi_HHHC,[r_C2H3; r_C2H4 ; r_C2H5; r_H3H4; r_H3H5; r_H4H5])+ sim(net_fphi_HHCBr,[r_C2H3; r_C2H4 ; r_C2Br6; r_H3H4; r_H3Br6; r_H4Br6])+...
%     sim(net_fphi_HHCBr,[r_C2H3; r_C2H5 ; r_C2Br6; r_H3H5; r_H3Br6; r_H5Br6])+...
%     sim(net_fphi_HHCBr,[r_C2H4; r_C2H5 ; r_C2Br6; r_H4H5; r_H4Br6; r_H5Br6])+...
%     sim(net_fphi_HHHBr,[r_H3H4; r_H3H5 ; r_H3Br6; r_H4H5; r_H4Br6; r_H5Br6]);


%%%%
ftheta_all= sim(net_ftheta_HHH,[r_H3H4; r_H3H5; r_H4H5]) + ...
    sim(net_ftheta_HHC,[r_C1H3; r_C1H4; r_H3H4]) + sim(net_ftheta_HHC,[r_C1H3; r_C1H5; r_H3H5]) + sim(net_ftheta_HHC,[r_C1H4; r_C1H5; r_H4H5])+...
    sim(net_ftheta_HHC,[r_C2H3; r_C2H4; r_H3H4]) + sim(net_ftheta_HHC,[r_C2H3; r_C2H5; r_H3H5]) + sim(net_ftheta_HHC,[r_C2H4; r_C2H5; r_H4H5])+...
    sim(net_ftheta_HHBr,[r_H3H4; r_H3Br6; r_H4Br6]) + sim(net_ftheta_HHBr,[r_H3H5; r_H3Br6; r_H5Br6]) + sim(net_ftheta_HHBr,[r_H4H5; r_H4Br6; r_H5Br6])+...
    sim(net_ftheta_HCC,[r_C1C2; r_C1H3; r_C2H3]) + sim(net_ftheta_HCC,[r_C1C2; r_C1H4; r_C2H4]) + sim(net_ftheta_HCC,[r_C1C2; r_C1H5; r_C2H5])+...
    sim(net_ftheta_HCBr,[r_C1H3; r_C1Br6; r_H3Br6]) + sim(net_ftheta_HCBr,[r_C1H4; r_C1Br6; r_H4Br6]) + sim(net_ftheta_HCBr,[r_C1H5; r_C1Br6; r_H5Br6])+...
    sim(net_ftheta_HCBr,[r_C2H3; r_C2Br6; r_H3Br6]) + sim(net_ftheta_HCBr,[r_C2H4; r_C2Br6; r_H4Br6]) + sim(net_ftheta_HCBr,[r_C2H5; r_C2Br6; r_H5Br6])+...
    sim(net_ftheta_CCBr,[r_C1C2; r_C1Br6; r_C2Br6]);



fphi_all= sim(net_fphi_HHHC,[r_C1H3; r_C1H4 ; r_C1H5; r_H3H4; r_H3H5; r_H4H5]) + sim(net_fphi_HHHC,[r_C2H3; r_C2H4 ; r_C2H5; r_H3H4; r_H3H5; r_H4H5]) + ...
    sim(net_fphi_HHHBr,[r_H3H4; r_H3H5 ; r_H3Br6; r_H4H5; r_H4Br6; r_H5Br6]) +...
    sim(net_fphi_HHCC,[r_C1C2; r_C1H3 ; r_C1H4; r_C2H3; r_C2H4; r_H3H4]) + sim(net_fphi_HHCC,[r_C1C2; r_C1H3 ; r_C1H5; r_C2H3; r_C2H5; r_H3H5]) + sim(net_fphi_HHCC,[r_C1C2; r_C1H4 ; r_C1H5; r_C2H4; r_C2H5; r_H4H5]) + ...
    sim(net_fphi_HHCBr,[r_C1H3; r_C1H4 ; r_C1Br6; r_H3H4; r_H3Br6; r_H4Br6]) + sim(net_fphi_HHCBr,[r_C1H3; r_C1H5 ; r_C1Br6; r_H3H5; r_H3Br6; r_H5Br6]) + sim(net_fphi_HHCBr,[r_C1H4; r_C1H5 ; r_C1Br6; r_H4H5; r_H4Br6; r_H5Br6]) + ...
    sim(net_fphi_HHCBr,[r_C2H3; r_C2H4 ; r_C2Br6; r_H3H4; r_H3Br6; r_H4Br6]) + sim(net_fphi_HHCBr,[r_C2H3; r_C2H5 ; r_C2Br6; r_H3H5; r_H3Br6; r_H5Br6]) + sim(net_fphi_HHCBr,[r_C2H4; r_C2H5 ; r_C2Br6; r_H4H5; r_H4Br6; r_H5Br6]) +...
    sim(net_fphi_HCCBr,[r_C1C2; r_C1H3 ; r_C1Br6; r_C2H3; r_C2Br6; r_H3Br6]) + sim(net_fphi_HCCBr,[r_C1C2; r_C1H4 ; r_C1Br6; r_C2H4; r_C2Br6; r_H4Br6]) +...
    sim(net_fphi_HCCBr,[r_C1C2; r_C1H5 ; r_C1Br6; r_C2H5; r_C2Br6; r_H5Br6]);

%%%%



Vhat_Milind= fr_all + ftheta_all + fphi_all;

Vhat_Milind = Vhat_Milind';


Ex_Milind=V-Vhat_Milind;
perf= sum(Ex_Milind.^2);



% Vhat_iQ= fc*(frHH + frHBr + fthetaHHBr);
% 
% Vhat(iQ)=Vhat_iQ; %#ok<AGROW>
% Ex(iQ,1)=V(iQ)-Vhat(iQ);  %#ok<AGROW>
% %SSE
% perf=perf+Ex(iQ,1)*Ex(iQ,1);

%%%%%$$$$$$$$$$@@@@@
                            
                            
% perf=0.0;
% 
% for iQ=1:1:Q
%     
%     %converting rs to my format
%     for i=1:1:clust_size(iQ)
%          for j=1:1:clust_size(iQ)
%              if i~=j
%                  r(i,j)=rs(i,j,iQ); %#ok<AGROW>
%              end
%          end
%     end
%     
%     
%     %calculating potential energy for iQ configuration
%      
%     frHH=0.0; frHBr=0.0; 
%     fthetaHHBr=0.0; 
%     
%     fc=1;%sim(net_fc,rij);
%     for i=1:1:clust_size(iQ)
%         for j=1:1:clust_size(iQ)
%             if j>i
%                 rij=r(i,j);
% %                 identify=horzcat(type(i),type(j));
%                 identify2=sprintf('%d %d', type(i),type(j));
%                 switch identify2
%                    case {'1 1'}
%                       frHH=frHH+sim(net_fr_HH,rij);
%                    case {'1 35', '35 1'}
%                       frHBr=frHBr+sim(net_fr_HBr,rij);
%                    otherwise
%                         disp(identify2)
%                    
%                 end
%                 
%             end
%             for k=1:1:clust_size(iQ)
%                 if (j>i && k>i && k>j)
%                     rik=r(i,k);
%                     rjk=r(j,k);
%                     identify3=sprintf('%d %d %d', type(i),type(j),type(k));
%                     switch identify3
%                        case {'1 1 35','1 35 1','35 1 1'}
%                           fthetaHHBr=fthetaHHBr+sim(net_ftheta_HHBr,[rij; rik ;rjk]);
%                        otherwise
%                         disp(identify3)
%                     end
%                                                       
%                 end
%                 
%             end
%            
%         end
%     end
%     
% %     
% %     for i=1:1:clust_size(iQ)
% %         for j=1:1:clust_size(iQ)
% %                                        
% %            
% %         end
% %     end
%     
% %     Vhat_iQ=fc*(frHH + frHBr + frHC + frCC + frCBr + ...
% %                 fthetaHHH + fthetaHHC + fthetaHHBr + fthetaHCBr + fthetaHCC + fthetaCCBr + ...
% %                 fphiHHHC + fphiHHHBr + fphiHHCC + fphiHHCBr + fphiHCCBr);
% Vhat_iQ=fc*(frHH + frHBr + fthetaHHBr);
%    
%     Vhat(iQ)=Vhat_iQ; %#ok<AGROW>
%     Ex(iQ,1)=V(iQ)-Vhat(iQ);  %#ok<AGROW>
%     %SSE
%     perf=perf+Ex(iQ,1)*Ex(iQ,1);
%      
%                  
%     
% end


i=1;%dummy line


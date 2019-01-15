

function [JE_Milind,JtJ_Milind,normJE_Milind]=calcjejj_GPES_Milind(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
                                net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
                                net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr,...
                                clust_size, Q, type, Ex,...
                     r_C1C2, r_C1H3, r_C1H4, r_C1H5, r_C1Br6, r_C2H3, r_C2H4, r_C2H5, r_C2Br6, r_H3H4, r_H3H5, r_H3Br6, r_H4H5, r_H4Br6, r_H5Br6)

% There will be 3 Jacobian matrices, one each for fc, fr, and ftheta

%%%% Jacobian matrix for fr
%J_fr will have dimensions of Q x ( 3N + 1 )
% First layer W (N x 1), b (N x 1), Second Layer W (1 x N ), b (1x1)

% W1HH = net_fr_HH.IW{1,1};
% W2HH = net_fr_HH.LW{2,1};
% B1HH = net_fr_HH.b{1};
% B2HH = net_fr_HH.b{2};
% 
% W1HBr = net_fr_HBr.IW{1,1};
% W2HBr = net_fr_HBr.LW{2,1};
% B1HBr = net_fr_HBr.b{1};
% B2HBr = net_fr_HBr.b{2};
% 
% S1 = 40;
% S2 = 1;


%---------------------
W1HH = net_fr_HH.IW{1,1};
W2HH = net_fr_HH.LW{2,1};
B1HH = net_fr_HH.b{1};
B2HH = net_fr_HH.b{2};

W1HBr = net_fr_HBr.IW{1,1};
W2HBr = net_fr_HBr.LW{2,1};
B1HBr = net_fr_HBr.b{1};
B2HBr = net_fr_HBr.b{2};

W1HC = net_fr_HC.IW{1,1};
W2HC = net_fr_HC.LW{2,1};
B1HC = net_fr_HC.b{1};
B2HC = net_fr_HC.b{2};

W1CC = net_fr_CC.IW{1,1};
W2CC = net_fr_CC.LW{2,1};
B1CC = net_fr_CC.b{1};
B2CC = net_fr_CC.b{2};

W1CBr = net_fr_CBr.IW{1,1};
W2CBr = net_fr_CBr.LW{2,1};
B1CBr = net_fr_CBr.b{1};
B2CBr = net_fr_CBr.b{2};


S1= 16;
S2= 1;

%%%%%$$$$$$$$$$@@@@@
% A1_H1H2 = logsig(W1HH*r_H1H2 + (repmat(B1HH*ones(1,(1)*1),1,Q)));
% A2_H1H2 = W2HH*A1_H1H2 + B2HH*ones(1,(1)*1);
% A1_H1H2 = kron(A1_H1H2,ones(1,S2));
% D2_H1H2 = nnmdlin(A2_H1H2);
% D1_H1H2 = nnmdlog(A1_H1H2,D2_H1H2,W2HH);
% jac1 = nnlmarq(kron(r_H1H2,ones(1,S2)),D1_H1H2);
% jac2 = nnlmarq(A1_H1H2,D2_H1H2);
% 
% jac_fr_H1H2 = [jac1,D1_H1H2',jac2,D2_H1H2']; 
% 
% 
% 
% 
% r_H1Br_H2Br= [r_H1Br; r_H2Br];
% 
% A1_H1Br = logsig(W1HBr*r_H1Br + (repmat(B1HBr*ones(1,(1)*1),1,Q)));
% A2_H1Br = W2HBr*A1_H1Br + B2HBr*ones(1,(1)*1);
% A1_H1Br = kron(A1_H1Br,ones(1,S2));
% D2_H1Br = nnmdlin(A2_H1Br);
% D1_H1Br = nnmdlog(A1_H1Br,D2_H1Br,W2HBr);
% jac1 = nnlmarq(kron(r_H1Br,ones(1,S2)),D1_H1Br);
% jac2 = nnlmarq(A1_H1Br,D2_H1Br);
% 
% jac_fr_H1Br = [jac1,D1_H1Br',jac2,D2_H1Br']; 
% 
% 
% 
% 
% A1_H2Br = logsig(W1HBr*r_H2Br + (repmat(B1HBr*ones(1,(1)*1),1,Q)));
% A2_H2Br = W2HBr*A1_H2Br + B2HBr*ones(1,(1)*1);
% A1_H2Br = kron(A1_H2Br,ones(1,S2));
% D2_H2Br = nnmdlin(A2_H2Br);
% D1_H2Br = nnmdlog(A1_H2Br,D2_H2Br,W2HBr);
% jac1 = nnlmarq(kron(r_H2Br,ones(1,S2)),D1_H2Br);
% jac2 = nnlmarq(A1_H2Br,D2_H2Br);
% 
% jac_fr_H2Br = [jac1,D1_H2Br',jac2,D2_H2Br']; 
% 
% 
% jac_fr_H1Br_H2Br= jac_fr_H1Br+jac_fr_H2Br;
% 
% J_fr_HH_HBr= [jac_fr_H1H2 jac_fr_H1Br_H2Br];

%-------------------------------------------------------
%-------------------------------------------------------
A1_C1C2 = logsig(W1CC*r_C1C2 + (repmat(B1CC*ones(1,(1)*1),1,Q)));
A2_C1C2 = W2CC*A1_C1C2 + B2CC*ones(1,(1)*1);
A1_C1C2 = kron(A1_C1C2,ones(1,S2));
D2_C1C2 = nnmdlin(A2_C1C2);
D1_C1C2 = nnmdlog(A1_C1C2,D2_C1C2,W2CC);
jac1 = nnlmarq(kron(r_C1C2,ones(1,S2)),D1_C1C2);
jac2 = nnlmarq(A1_C1C2,D2_C1C2);

jac_fr_C1C2 = [jac1,D1_C1C2',jac2,D2_C1C2']; 



A1_C1H3 = logsig(W1HC*r_C1H3 + (repmat(B1HC*ones(1,(1)*1),1,Q)));
A2_C1H3 = W2HC*A1_C1H3 + B2HC*ones(1,(1)*1);
A1_C1H3 = kron(A1_C1H3,ones(1,S2));
D2_C1H3 = nnmdlin(A2_C1H3);
D1_C1H3 = nnmdlog(A1_C1H3,D2_C1H3,W2HC);
jac1 = nnlmarq(kron(r_C1H3,ones(1,S2)),D1_C1H3);
jac2 = nnlmarq(A1_C1H3,D2_C1H3);

jac_fr_C1H3 = [jac1,D1_C1H3',jac2,D2_C1H3']; 



A1_C1H4 = logsig(W1HC*r_C1H4 + (repmat(B1HC*ones(1,(1)*1),1,Q)));
A2_C1H4 = W2HC*A1_C1H4 + B2HC*ones(1,(1)*1);
A1_C1H4 = kron(A1_C1H4,ones(1,S2));
D2_C1H4 = nnmdlin(A2_C1H4);
D1_C1H4 = nnmdlog(A1_C1H4,D2_C1H4,W2HC);
jac1 = nnlmarq(kron(r_C1H4,ones(1,S2)),D1_C1H4);
jac2 = nnlmarq(A1_C1H4,D2_C1H4);

jac_fr_C1H4 = [jac1,D1_C1H4',jac2,D2_C1H4']; 



A1_C1H5 = logsig(W1HC*r_C1H5 + (repmat(B1HC*ones(1,(1)*1),1,Q)));
A2_C1H5 = W2HC*A1_C1H5 + B2HC*ones(1,(1)*1);
A1_C1H5 = kron(A1_C1H5,ones(1,S2));
D2_C1H5 = nnmdlin(A2_C1H5);
D1_C1H5 = nnmdlog(A1_C1H5,D2_C1H5,W2HC);
jac1 = nnlmarq(kron(r_C1H5,ones(1,S2)),D1_C1H5);
jac2 = nnlmarq(A1_C1H5,D2_C1H5);

jac_fr_C1H5 = [jac1,D1_C1H5',jac2,D2_C1H5']; 



A1_C1Br6 = logsig(W1CBr*r_C1Br6 + (repmat(B1CBr*ones(1,(1)*1),1,Q)));
A2_C1Br6 = W2CBr*A1_C1Br6 + B2CBr*ones(1,(1)*1);
A1_C1Br6 = kron(A1_C1Br6,ones(1,S2));
D2_C1Br6 = nnmdlin(A2_C1Br6);
D1_C1Br6 = nnmdlog(A1_C1Br6,D2_C1Br6,W2CBr);
jac1 = nnlmarq(kron(r_C1Br6,ones(1,S2)),D1_C1Br6);
jac2 = nnlmarq(A1_C1Br6,D2_C1Br6);

jac_fr_C1Br6 = [jac1,D1_C1Br6',jac2,D2_C1Br6']; 



A1_C2H3 = logsig(W1HC*r_C2H3 + (repmat(B1HC*ones(1,(1)*1),1,Q)));
A2_C2H3 = W2HC*A1_C2H3 + B2HC*ones(1,(1)*1);
A1_C2H3 = kron(A1_C2H3,ones(1,S2));
D2_C2H3 = nnmdlin(A2_C2H3);
D1_C2H3 = nnmdlog(A1_C2H3,D2_C2H3,W2HC);
jac1 = nnlmarq(kron(r_C2H3,ones(1,S2)),D1_C2H3);
jac2 = nnlmarq(A1_C2H3,D2_C2H3);

jac_fr_C2H3 = [jac1,D1_C2H3',jac2,D2_C2H3']; 




A1_C2H4 = logsig(W1HC*r_C2H4 + (repmat(B1HC*ones(1,(1)*1),1,Q)));
A2_C2H4 = W2HC*A1_C2H4 + B2HC*ones(1,(1)*1);
A1_C2H4 = kron(A1_C2H4,ones(1,S2));
D2_C2H4 = nnmdlin(A2_C2H4);
D1_C2H4 = nnmdlog(A1_C2H4,D2_C2H4,W2HC);
jac1 = nnlmarq(kron(r_C2H4,ones(1,S2)),D1_C2H4);
jac2 = nnlmarq(A1_C2H4,D2_C2H4);

jac_fr_C2H4 = [jac1,D1_C2H4',jac2,D2_C2H4']; 




A1_C2H5 = logsig(W1HC*r_C2H5 + (repmat(B1HC*ones(1,(1)*1),1,Q)));
A2_C2H5 = W2HC*A1_C2H5 + B2HC*ones(1,(1)*1);
A1_C2H5 = kron(A1_C2H5,ones(1,S2));
D2_C2H5 = nnmdlin(A2_C2H5);
D1_C2H5 = nnmdlog(A1_C2H5,D2_C2H5,W2HC);
jac1 = nnlmarq(kron(r_C2H5,ones(1,S2)),D1_C2H5);
jac2 = nnlmarq(A1_C2H5,D2_C2H5);

jac_fr_C2H5 = [jac1,D1_C2H5',jac2,D2_C2H5']; 



A1_C2Br6 = logsig(W1CBr*r_C2Br6 + (repmat(B1CBr*ones(1,(1)*1),1,Q)));
A2_C2Br6 = W2CBr*A1_C2Br6 + B2CBr*ones(1,(1)*1);
A1_C2Br6 = kron(A1_C2Br6,ones(1,S2));
D2_C2Br6 = nnmdlin(A2_C2Br6);
D1_C2Br6 = nnmdlog(A1_C2Br6,D2_C2Br6,W2CBr);
jac1 = nnlmarq(kron(r_C2Br6,ones(1,S2)),D1_C2Br6);
jac2 = nnlmarq(A1_C2Br6,D2_C2Br6);

jac_fr_C2Br6 = [jac1,D1_C2Br6',jac2,D2_C2Br6']; 




A1_H3H4 = logsig(W1HH*r_H3H4 + (repmat(B1HH*ones(1,(1)*1),1,Q)));
A2_H3H4 = W2HH*A1_H3H4 + B2HH*ones(1,(1)*1);
A1_H3H4 = kron(A1_H3H4,ones(1,S2));
D2_H3H4 = nnmdlin(A2_H3H4);
D1_H3H4 = nnmdlog(A1_H3H4,D2_H3H4,W2HH);
jac1 = nnlmarq(kron(r_H3H4,ones(1,S2)),D1_H3H4);
jac2 = nnlmarq(A1_H3H4,D2_H3H4);

jac_fr_H3H4 = [jac1,D1_H3H4',jac2,D2_H3H4']; 




A1_H3H5 = logsig(W1HH*r_H3H5 + (repmat(B1HH*ones(1,(1)*1),1,Q)));
A2_H3H5 = W2HH*A1_H3H5 + B2HH*ones(1,(1)*1);
A1_H3H5 = kron(A1_H3H5,ones(1,S2));
D2_H3H5 = nnmdlin(A2_H3H5);
D1_H3H5 = nnmdlog(A1_H3H5,D2_H3H5,W2HH);
jac1 = nnlmarq(kron(r_H3H5,ones(1,S2)),D1_H3H5);
jac2 = nnlmarq(A1_H3H5,D2_H3H5);

jac_fr_H3H5 = [jac1,D1_H3H5',jac2,D2_H3H5']; 




A1_H3Br6 = logsig(W1HBr*r_H3Br6 + (repmat(B1HBr*ones(1,(1)*1),1,Q)));
A2_H3Br6 = W2HBr*A1_H3Br6 + B2HBr*ones(1,(1)*1);
A1_H3Br6 = kron(A1_H3Br6,ones(1,S2));
D2_H3Br6 = nnmdlin(A2_H3Br6);
D1_H3Br6 = nnmdlog(A1_H3Br6,D2_H3Br6,W2HBr);
jac1 = nnlmarq(kron(r_H3Br6,ones(1,S2)),D1_H3Br6);
jac2 = nnlmarq(A1_H3Br6,D2_H3Br6);

jac_fr_H3Br6 = [jac1,D1_H3Br6',jac2,D2_H3Br6']; 




A1_H4H5 = logsig(W1HH*r_H4H5 + (repmat(B1HH*ones(1,(1)*1),1,Q)));
A2_H4H5 = W2HH*A1_H4H5 + B2HH*ones(1,(1)*1);
A1_H4H5 = kron(A1_H4H5,ones(1,S2));
D2_H4H5 = nnmdlin(A2_H4H5);
D1_H4H5 = nnmdlog(A1_H4H5,D2_H4H5,W2HH);
jac1 = nnlmarq(kron(r_H4H5,ones(1,S2)),D1_H4H5);
jac2 = nnlmarq(A1_H4H5,D2_H4H5);

jac_fr_H4H5 = [jac1,D1_H4H5',jac2,D2_H4H5']; 



A1_H4Br6 = logsig(W1HBr*r_H4Br6 + (repmat(B1HBr*ones(1,(1)*1),1,Q)));
A2_H4Br6 = W2HBr*A1_H4Br6 + B2HBr*ones(1,(1)*1);
A1_H4Br6 = kron(A1_H4Br6,ones(1,S2));
D2_H4Br6 = nnmdlin(A2_H4Br6);
D1_H4Br6 = nnmdlog(A1_H4Br6,D2_H4Br6,W2HBr);
jac1 = nnlmarq(kron(r_H4Br6,ones(1,S2)),D1_H4Br6);
jac2 = nnlmarq(A1_H4Br6,D2_H4Br6);

jac_fr_H4Br6 = [jac1,D1_H4Br6',jac2,D2_H4Br6']; 



A1_H5Br6 = logsig(W1HBr*r_H5Br6 + (repmat(B1HBr*ones(1,(1)*1),1,Q)));
A2_H5Br6 = W2HBr*A1_H5Br6 + B2HBr*ones(1,(1)*1);
A1_H5Br6 = kron(A1_H5Br6,ones(1,S2));
D2_H5Br6 = nnmdlin(A2_H5Br6);
D1_H5Br6 = nnmdlog(A1_H5Br6,D2_H5Br6,W2HBr);
jac1 = nnlmarq(kron(r_H5Br6,ones(1,S2)),D1_H5Br6);
jac2 = nnlmarq(A1_H5Br6,D2_H5Br6);

jac_fr_H5Br6 = [jac1,D1_H5Br6',jac2,D2_H5Br6']; 




jac_fr_HH= jac_fr_H3H4 + jac_fr_H3H5 + jac_fr_H4H5;
jac_fr_HC= jac_fr_C1H3 + jac_fr_C1H4 + jac_fr_C1H5 + jac_fr_C2H3 + jac_fr_C2H4 + jac_fr_C2H5;
jac_fr_HBr= jac_fr_H3Br6 + jac_fr_H4Br6 + jac_fr_H5Br6;
jac_fr_CC= jac_fr_C1C2;
jac_fr_CBr= jac_fr_C1Br6 + jac_fr_C2Br6;

J_fr_all= [jac_fr_HH jac_fr_HC jac_fr_HBr jac_fr_CC jac_fr_CBr];

clear jac_fr_HH jac_fr_HC jac_fr_HBr jac_fr_CC jac_fr_CBr;
clear W1HH W2HH B1HH B2HH W1HBr W2HBr B1HBr B2HBr W1HC W2HC B1HC B2HC W1CC W2CC B1CC B2CC W1CBr W2CBr B1CBr B2CBr 

clear A1_C1C2 A2_C1C2 A1_C1C2 D2_C1C2 D1_C1C2 jac1 jac2 jac_fr_C1C2 	
clear A1_C1H3 A2_C1H3 A1_C1H3 D2_C1H3 D1_C1H3 jac1 jac2 jac_fr_C1H3 				
clear A1_C1H4 A2_C1H4 A1_C1H4 D2_C1H4 D1_C1H4 jac1 jac2 jac_fr_C1H4 				
clear A1_C1H5 	A2_C1H5 	A1_C1H5 	D2_C1H5 	D1_C1H5 	jac1 	jac2 		jac_fr_C1H5 				
clear A1_C1Br6 	A2_C1Br6 	A1_C1Br6 	D2_C1Br6 	D1_C1Br6 	jac1 	jac2 		jac_fr_C1Br6
clear A1_C2H3 	A2_C2H3 	A1_C2H3 	D2_C2H3 	D1_C2H3 	jac1 	jac2 		jac_fr_C2H3 
clear A1_C2H4 	A2_C2H4 	A1_C2H4 	D2_C2H4 	D1_C2H4 	jac1 	jac2 		jac_fr_C2H4 
clear A1_C2H5 	A2_C2H5 	A1_C2H5 	D2_C2H5 	D1_C2H5 	jac1 	jac2 		jac_fr_C2H5
clear A1_C2Br6 	A2_C2Br6 	A1_C2Br6 	D2_C2Br6 	D1_C2Br6 	jac1 	jac2 		jac_fr_C2Br6 
clear A1_H3H4 	A2_H3H4 	A1_H3H4 	D2_H3H4 	D1_H3H4 	jac1 	jac2 		jac_fr_H3H4 
clear A1_H3H5 	A2_H3H5 	A1_H3H5 	D2_H3H5 	D1_H3H5 	jac1 	jac2 		jac_fr_H3H5 	
clear A1_H3Br6 	A2_H3Br6 	A1_H3Br6 	D2_H3Br6 	D1_H3Br6 	jac1 	jac2 		jac_fr_H3Br6 	
clear A1_H4H5 	A2_H4H5 	A1_H4H5 	D2_H4H5 	D1_H4H5 	jac1 	jac2 		jac_fr_H4H5 
clear A1_H4Br6 	A2_H4Br6 	A1_H4Br6 	D2_H4Br6 	D1_H4Br6 	jac1 	jac2 		jac_fr_H4Br6 
clear A1_H5Br6 	A2_H5Br6 	A1_H5Br6 	D2_H5Br6 	D1_H5Br6 	jac1 	jac2 		jac_fr_H5Br6 


%-------------------------------------------------------
%-------------------------------------------------------



 

W1HHH = net_ftheta_HHH.IW{1,1};
W2HHH = net_ftheta_HHH.LW{2,1};
B1HHH = net_ftheta_HHH.b{1};
B2HHH = net_ftheta_HHH.b{2};

W1HHC = net_ftheta_HHC.IW{1,1};
W2HHC = net_ftheta_HHC.LW{2,1};
B1HHC = net_ftheta_HHC.b{1};
B2HHC = net_ftheta_HHC.b{2};

W1HHBr = net_ftheta_HHBr.IW{1,1};
W2HHBr = net_ftheta_HHBr.LW{2,1};
B1HHBr = net_ftheta_HHBr.b{1};
B2HHBr = net_ftheta_HHBr.b{2};

W1HCBr = net_ftheta_HCBr.IW{1,1};
W2HCBr = net_ftheta_HCBr.LW{2,1};
B1HCBr = net_ftheta_HCBr.b{1};
B2HCBr = net_ftheta_HCBr.b{2};

W1HCC = net_ftheta_HCC.IW{1,1};
W2HCC = net_ftheta_HCC.LW{2,1};
B1HCC = net_ftheta_HCC.b{1};
B2HCC = net_ftheta_HCC.b{2};

W1CCBr = net_ftheta_CCBr.IW{1,1};
W2CCBr = net_ftheta_CCBr.LW{2,1};
B1CCBr = net_ftheta_CCBr.b{1};
B2CCBr = net_ftheta_CCBr.b{2};



% W1HHBr = net_ftheta_HHBr.IW{1,1};
% W2HHBr = net_ftheta_HHBr.LW{2,1};
% B1HHBr = net_ftheta_HHBr.b{1};
% B2HHBr = net_ftheta_HHBr.b{2};

S1 = 24;
S2 = 1;


A1_C1C2H3 = logsig(W1HCC*[r_C1C2; r_C1H3; r_C2H3] + (repmat(B1HCC*ones(1,(1)*1),1,Q ) ));
A2_C1C2H3 = W2HCC*A1_C1C2H3 + B2HCC*ones(1,(1)*1);
A1_C1C2H3 = kron(A1_C1C2H3,ones(1,S2));
D2_C1C2H3 = nnmdlin(A2_C1C2H3);
D1_C1C2H3 = nnmdlog(A1_C1C2H3,D2_C1C2H3,W2HCC);
jac1 = nnlmarq(kron([r_C1C2; r_C1H3; r_C2H3],ones(1,S2)),D1_C1C2H3);
jac2 = nnlmarq(A1_C1C2H3,D2_C1C2H3);

jac_ftheta_C1C2H3 = [jac1,D1_C1C2H3',jac2,D2_C1C2H3']; 



A1_C1C2H4 = logsig(W1HCC*[r_C1C2; r_C1H4; r_C2H4] + (repmat(B1HCC*ones(1,(1)*1),1,Q ) ));
A2_C1C2H4 = W2HCC*A1_C1C2H4 + B2HCC*ones(1,(1)*1);
A1_C1C2H4 = kron(A1_C1C2H4,ones(1,S2));
D2_C1C2H4 = nnmdlin(A2_C1C2H4);
D1_C1C2H4 = nnmdlog(A1_C1C2H4,D2_C1C2H4,W2HCC);
jac1 = nnlmarq(kron([r_C1C2; r_C1H4; r_C2H4],ones(1,S2)),D1_C1C2H4);
jac2 = nnlmarq(A1_C1C2H4,D2_C1C2H4);

jac_ftheta_C1C2H4 = [jac1,D1_C1C2H4',jac2,D2_C1C2H4']; 



A1_C1C2H5 = logsig(W1HCC*[r_C1C2; r_C1H5; r_C2H5] + (repmat(B1HCC*ones(1,(1)*1),1,Q ) ));
A2_C1C2H5 = W2HCC*A1_C1C2H5 + B2HCC*ones(1,(1)*1);
A1_C1C2H5 = kron(A1_C1C2H5,ones(1,S2));
D2_C1C2H5 = nnmdlin(A2_C1C2H5);
D1_C1C2H5 = nnmdlog(A1_C1C2H5,D2_C1C2H5,W2HCC);
jac1 = nnlmarq(kron([r_C1C2; r_C1H5; r_C2H5],ones(1,S2)),D1_C1C2H5);
jac2 = nnlmarq(A1_C1C2H5,D2_C1C2H5);

jac_ftheta_C1C2H5 = [jac1,D1_C1C2H5',jac2,D2_C1C2H5']; 



A1_C1C2Br6 = logsig(W1CCBr*[r_C1C2; r_C1Br6; r_C2Br6] + (repmat(B1CCBr*ones(1,(1)*1),1,Q ) ));
A2_C1C2Br6 = W2CCBr*A1_C1C2Br6 + B2CCBr*ones(1,(1)*1);
A1_C1C2Br6 = kron(A1_C1C2Br6,ones(1,S2));
D2_C1C2Br6 = nnmdlin(A2_C1C2Br6);
D1_C1C2Br6 = nnmdlog(A1_C1C2Br6,D2_C1C2Br6,W2CCBr);
jac1 = nnlmarq(kron([r_C1C2; r_C1Br6; r_C2Br6],ones(1,S2)),D1_C1C2Br6);
jac2 = nnlmarq(A1_C1C2Br6,D2_C1C2Br6);

jac_ftheta_C1C2Br6 = [jac1,D1_C1C2Br6',jac2,D2_C1C2Br6']; 

%%%%---

A1_C1H3H4 = logsig(W1HHC*[r_C1H3; r_C1H4; r_H3H4] + (repmat(B1HHC*ones(1,(1)*1),1,Q ) ));
A2_C1H3H4 = W2HHC*A1_C1H3H4 + B2HHC*ones(1,(1)*1);
A1_C1H3H4 = kron(A1_C1H3H4,ones(1,S2));
D2_C1H3H4 = nnmdlin(A2_C1H3H4);
D1_C1H3H4 = nnmdlog(A1_C1H3H4,D2_C1H3H4,W2HHC);
jac1 = nnlmarq(kron([r_C1H3; r_C1H4; r_H3H4],ones(1,S2)),D1_C1H3H4);
jac2 = nnlmarq(A1_C1H3H4,D2_C1H3H4);

jac_ftheta_C1H3H4 = [jac1,D1_C1H3H4',jac2,D2_C1H3H4']; 



A1_C1H3H5 = logsig(W1HHC*[r_C1H3; r_C1H5; r_H3H5] + (repmat(B1HHC*ones(1,(1)*1),1,Q ) ));
A2_C1H3H5 = W2HHC*A1_C1H3H5 + B2HHC*ones(1,(1)*1);
A1_C1H3H5 = kron(A1_C1H3H5,ones(1,S2));
D2_C1H3H5 = nnmdlin(A2_C1H3H5);
D1_C1H3H5 = nnmdlog(A1_C1H3H5,D2_C1H3H5,W2HHC);
jac1 = nnlmarq(kron([r_C1H3; r_C1H5; r_H3H5],ones(1,S2)),D1_C1H3H5);
jac2 = nnlmarq(A1_C1H3H5,D2_C1H3H5);

jac_ftheta_C1H3H5 = [jac1,D1_C1H3H5',jac2,D2_C1H3H5']; 



A1_C1H4H5 = logsig(W1HHC*[r_C1H4; r_C1H5; r_H4H5] + (repmat(B1HHC*ones(1,(1)*1),1,Q ) ));
A2_C1H4H5 = W2HHC*A1_C1H4H5 + B2HHC*ones(1,(1)*1);
A1_C1H4H5 = kron(A1_C1H4H5,ones(1,S2));
D2_C1H4H5 = nnmdlin(A2_C1H4H5);
D1_C1H4H5 = nnmdlog(A1_C1H4H5,D2_C1H4H5,W2HHC);
jac1 = nnlmarq(kron([r_C1H4; r_C1H5; r_H4H5],ones(1,S2)),D1_C1H4H5);
jac2 = nnlmarq(A1_C1H4H5,D2_C1H4H5);

jac_ftheta_C1H4H5 = [jac1,D1_C1H4H5',jac2,D2_C1H4H5']; 





A1_C2H3H4 = logsig(W1HHC*[r_C2H3; r_C2H4; r_H3H4] + (repmat(B1HHC*ones(1,(1)*1),1,Q ) ));
A2_C2H3H4 = W2HHC*A1_C2H3H4 + B2HHC*ones(1,(1)*1);
A1_C2H3H4 = kron(A1_C2H3H4,ones(1,S2));
D2_C2H3H4 = nnmdlin(A2_C2H3H4);
D1_C2H3H4 = nnmdlog(A1_C2H3H4,D2_C2H3H4,W2HHC);
jac1 = nnlmarq(kron([r_C2H3; r_C2H4; r_H3H4],ones(1,S2)),D1_C2H3H4);
jac2 = nnlmarq(A1_C2H3H4,D2_C2H3H4);

jac_ftheta_C2H3H4 = [jac1,D1_C2H3H4',jac2,D2_C2H3H4']; 



A1_C2H3H5 = logsig(W1HHC*[r_C2H3; r_C2H5; r_H3H5] + (repmat(B1HHC*ones(1,(1)*1),1,Q ) ));
A2_C2H3H5 = W2HHC*A1_C2H3H5 + B2HHC*ones(1,(1)*1);
A1_C2H3H5 = kron(A1_C2H3H5,ones(1,S2));
D2_C2H3H5 = nnmdlin(A2_C2H3H5);
D1_C2H3H5 = nnmdlog(A1_C2H3H5,D2_C2H3H5,W2HHC);
jac1 = nnlmarq(kron([r_C2H3; r_C2H5; r_H3H5],ones(1,S2)),D1_C2H3H5);
jac2 = nnlmarq(A1_C2H3H5,D2_C2H3H5);

jac_ftheta_C2H3H5 = [jac1,D1_C2H3H5',jac2,D2_C2H3H5']; 



A1_C2H4H5 = logsig(W1HHC*[r_C2H4; r_C2H5; r_H4H5] + (repmat(B1HHC*ones(1,(1)*1),1,Q ) ));
A2_C2H4H5 = W2HHC*A1_C2H4H5 + B2HHC*ones(1,(1)*1);
A1_C2H4H5 = kron(A1_C2H4H5,ones(1,S2));
D2_C2H4H5 = nnmdlin(A2_C2H4H5);
D1_C2H4H5 = nnmdlog(A1_C2H4H5,D2_C2H4H5,W2HHC);
jac1 = nnlmarq(kron([r_C2H4; r_C2H5; r_H4H5],ones(1,S2)),D1_C2H4H5);
jac2 = nnlmarq(A1_C2H4H5,D2_C2H4H5);

jac_ftheta_C2H4H5 = [jac1,D1_C2H4H5',jac2,D2_C2H4H5']; 




A1_C1H3Br6 = logsig(W1HCBr*[r_C1H3; r_C1Br6; r_H3Br6] + (repmat(B1HCBr*ones(1,(1)*1),1,Q ) ));
A2_C1H3Br6 = W2HCBr*A1_C1H3Br6 + B2HCBr*ones(1,(1)*1);
A1_C1H3Br6 = kron(A1_C1H3Br6,ones(1,S2));
D2_C1H3Br6 = nnmdlin(A2_C1H3Br6);
D1_C1H3Br6 = nnmdlog(A1_C1H3Br6,D2_C1H3Br6,W2HCBr);
jac1 = nnlmarq(kron([r_C1H3; r_C1Br6; r_H3Br6],ones(1,S2)),D1_C1H3Br6);
jac2 = nnlmarq(A1_C1H3Br6,D2_C1H3Br6);

jac_ftheta_C1H3Br6 = [jac1,D1_C1H3Br6',jac2,D2_C1H3Br6']; 



A1_C1H4Br6 = logsig(W1HCBr*[r_C1H4; r_C1Br6; r_H4Br6] + (repmat(B1HCBr*ones(1,(1)*1),1,Q ) ));
A2_C1H4Br6 = W2HCBr*A1_C1H4Br6 + B2HCBr*ones(1,(1)*1);
A1_C1H4Br6 = kron(A1_C1H4Br6,ones(1,S2));
D2_C1H4Br6 = nnmdlin(A2_C1H4Br6);
D1_C1H4Br6 = nnmdlog(A1_C1H4Br6,D2_C1H4Br6,W2HCBr);
jac1 = nnlmarq(kron([r_C1H4; r_C1Br6; r_H4Br6],ones(1,S2)),D1_C1H4Br6);
jac2 = nnlmarq(A1_C1H4Br6,D2_C1H4Br6);

jac_ftheta_C1H4Br6 = [jac1,D1_C1H4Br6',jac2,D2_C1H4Br6']; 




A1_C1H5Br6 = logsig(W1HCBr*[r_C1H5; r_C1Br6; r_H5Br6] + (repmat(B1HCBr*ones(1,(1)*1),1,Q ) ));
A2_C1H5Br6 = W2HCBr*A1_C1H5Br6 + B2HCBr*ones(1,(1)*1);
A1_C1H5Br6 = kron(A1_C1H5Br6,ones(1,S2));
D2_C1H5Br6 = nnmdlin(A2_C1H5Br6);
D1_C1H5Br6 = nnmdlog(A1_C1H5Br6,D2_C1H5Br6,W2HCBr);
jac1 = nnlmarq(kron([r_C1H5; r_C1Br6; r_H5Br6],ones(1,S2)),D1_C1H5Br6);
jac2 = nnlmarq(A1_C1H5Br6,D2_C1H5Br6);

jac_ftheta_C1H5Br6 = [jac1,D1_C1H5Br6',jac2,D2_C1H5Br6']; 





A1_C2H3Br6 = logsig(W1HCBr*[r_C2H3; r_C2Br6; r_H3Br6] + (repmat(B1HCBr*ones(1,(1)*1),1,Q ) ));
A2_C2H3Br6 = W2HCBr*A1_C2H3Br6 + B2HCBr*ones(1,(1)*1);
A1_C2H3Br6 = kron(A1_C2H3Br6,ones(1,S2));
D2_C2H3Br6 = nnmdlin(A2_C2H3Br6);
D1_C2H3Br6 = nnmdlog(A1_C2H3Br6,D2_C2H3Br6,W2HCBr);
jac1 = nnlmarq(kron([r_C2H3; r_C2Br6; r_H3Br6],ones(1,S2)),D1_C2H3Br6);
jac2 = nnlmarq(A1_C2H3Br6,D2_C2H3Br6);

jac_ftheta_C2H3Br6 = [jac1,D1_C2H3Br6',jac2,D2_C2H3Br6']; 



A1_C2H4Br6 = logsig(W1HCBr*[r_C2H4; r_C2Br6; r_H4Br6] + (repmat(B1HCBr*ones(1,(1)*1),1,Q ) ));
A2_C2H4Br6 = W2HCBr*A1_C2H4Br6 + B2HCBr*ones(1,(1)*1);
A1_C2H4Br6 = kron(A1_C2H4Br6,ones(1,S2));
D2_C2H4Br6 = nnmdlin(A2_C2H4Br6);
D1_C2H4Br6 = nnmdlog(A1_C2H4Br6,D2_C2H4Br6,W2HCBr);
jac1 = nnlmarq(kron([r_C2H4; r_C2Br6; r_H4Br6],ones(1,S2)),D1_C2H4Br6);
jac2 = nnlmarq(A1_C2H4Br6,D2_C2H4Br6);

jac_ftheta_C2H4Br6 = [jac1,D1_C2H4Br6',jac2,D2_C2H4Br6']; 




A1_C2H5Br6 = logsig(W1HCBr*[r_C2H5; r_C2Br6; r_H5Br6] + (repmat(B1HCBr*ones(1,(1)*1),1,Q ) ));
A2_C2H5Br6 = W2HCBr*A1_C2H5Br6 + B2HCBr*ones(1,(1)*1);
A1_C2H5Br6 = kron(A1_C2H5Br6,ones(1,S2));
D2_C2H5Br6 = nnmdlin(A2_C2H5Br6);
D1_C2H5Br6 = nnmdlog(A1_C2H5Br6,D2_C2H5Br6,W2HCBr);
jac1 = nnlmarq(kron([r_C2H5; r_C2Br6; r_H5Br6],ones(1,S2)),D1_C2H5Br6);
jac2 = nnlmarq(A1_C2H5Br6,D2_C2H5Br6);

jac_ftheta_C2H5Br6 = [jac1,D1_C2H5Br6',jac2,D2_C2H5Br6']; 




A1_H3H4H5 = logsig(W1HHH*[r_H3H4; r_H3H5; r_H4H5] + (repmat(B1HHH*ones(1,(1)*1),1,Q ) ));
A2_H3H4H5 = W2HHH*A1_H3H4H5 + B2HHH*ones(1,(1)*1);
A1_H3H4H5 = kron(A1_H3H4H5,ones(1,S2));
D2_H3H4H5 = nnmdlin(A2_H3H4H5);
D1_H3H4H5 = nnmdlog(A1_H3H4H5,D2_H3H4H5,W2HHH);
jac1 = nnlmarq(kron([r_H3H4; r_H3H5; r_H4H5],ones(1,S2)),D1_H3H4H5);
jac2 = nnlmarq(A1_H3H4H5,D2_H3H4H5);

jac_ftheta_H3H4H5 = [jac1,D1_H3H4H5',jac2,D2_H3H4H5']; 




A1_H3H4Br6 = logsig(W1HHBr*[r_H3H4; r_H3Br6; r_H4Br6] + (repmat(B1HHBr*ones(1,(1)*1),1,Q ) ));
A2_H3H4Br6 = W2HHBr*A1_H3H4Br6 + B2HHBr*ones(1,(1)*1);
A1_H3H4Br6 = kron(A1_H3H4Br6,ones(1,S2));
D2_H3H4Br6 = nnmdlin(A2_H3H4Br6);
D1_H3H4Br6 = nnmdlog(A1_H3H4Br6,D2_H3H4Br6,W2HHBr);
jac1 = nnlmarq(kron([r_H3H4; r_H3Br6; r_H4Br6],ones(1,S2)),D1_H3H4Br6);
jac2 = nnlmarq(A1_H3H4Br6,D2_H3H4Br6);

jac_ftheta_H3H4Br6 = [jac1,D1_H3H4Br6',jac2,D2_H3H4Br6']; 




A1_H3H5Br6 = logsig(W1HHBr*[r_H3H5; r_H3Br6; r_H5Br6] + (repmat(B1HHBr*ones(1,(1)*1),1,Q ) ));
A2_H3H5Br6 = W2HHBr*A1_H3H5Br6 + B2HHBr*ones(1,(1)*1);
A1_H3H5Br6 = kron(A1_H3H5Br6,ones(1,S2));
D2_H3H5Br6 = nnmdlin(A2_H3H5Br6);
D1_H3H5Br6 = nnmdlog(A1_H3H5Br6,D2_H3H5Br6,W2HHBr);
jac1 = nnlmarq(kron([r_H3H5; r_H3Br6; r_H5Br6],ones(1,S2)),D1_H3H5Br6);
jac2 = nnlmarq(A1_H3H5Br6,D2_H3H5Br6);

jac_ftheta_H3H5Br6 = [jac1,D1_H3H5Br6',jac2,D2_H3H5Br6']; 



A1_H4H5Br6 = logsig(W1HHBr*[r_H4H5; r_H4Br6; r_H5Br6] + (repmat(B1HHBr*ones(1,(1)*1),1,Q ) ));
A2_H4H5Br6 = W2HHBr*A1_H4H5Br6 + B2HHBr*ones(1,(1)*1);
A1_H4H5Br6 = kron(A1_H4H5Br6,ones(1,S2));
D2_H4H5Br6 = nnmdlin(A2_H4H5Br6);
D1_H4H5Br6 = nnmdlog(A1_H4H5Br6,D2_H4H5Br6,W2HHBr);
jac1 = nnlmarq(kron([r_H4H5; r_H4Br6; r_H5Br6],ones(1,S2)),D1_H4H5Br6);
jac2 = nnlmarq(A1_H4H5Br6,D2_H4H5Br6);

jac_ftheta_H4H5Br6 = [jac1,D1_H4H5Br6',jac2,D2_H4H5Br6']; 





jac_ftheta_HHH= [jac_ftheta_H3H4H5];
jac_ftheta_HHC= [jac_ftheta_C1H3H4 + jac_ftheta_C1H3H5 + jac_ftheta_C1H4H5 + jac_ftheta_C2H3H4 + jac_ftheta_C2H3H5 + jac_ftheta_C2H4H5];
jac_ftheta_HHBr= [jac_ftheta_H3H4Br6 + jac_ftheta_H3H5Br6 + jac_ftheta_H4H5Br6];
jac_ftheta_HCC= [jac_ftheta_C1C2H3 + jac_ftheta_C1C2H4 + jac_ftheta_C1C2H5];
jac_ftheta_HCBr= [jac_ftheta_C1H3Br6 + jac_ftheta_C1H4Br6 + jac_ftheta_C1H5Br6 + jac_ftheta_C2H3Br6 + jac_ftheta_C2H4Br6 + jac_ftheta_C2H5Br6];
jac_ftheta_CCBr= [jac_ftheta_C1C2Br6];


J_ftheta_all= [jac_ftheta_HHH jac_ftheta_HHC jac_ftheta_HHBr jac_ftheta_HCC jac_ftheta_HCBr jac_ftheta_CCBr];

clear jac_ftheta_HHH jac_ftheta_HHC jac_ftheta_HHBr jac_ftheta_HCC jac_ftheta_HCBr jac_ftheta_CCBr;
clear W1HHH 	W2HHH 	B1HHH 	B2HHH 		W1HHC 	W2HHC 	B1HHC 	B2HHC 		W1HHBr 	W2HHBr 	B1HHBr 	B2HHBr 		W1HCBr 	W2HCBr 	B1HCBr 	B2HCBr 		W1HCC 	W2HCC 	B1HCC 	B2HCC 		W1CCBr 	W2CCBr 	B1CCBr 	B2CCBr 

clear A1_C1C2H3 	A2_C1C2H3 	A1_C1C2H3 	D2_C1C2H3 	D1_C1C2H3 	jac1 	jac2 		jac_ftheta_C1C2H3 				A1_C1C2H4 	A2_C1C2H4 	A1_C1C2H4 	D2_C1C2H4 	D1_C1C2H4 	jac1 	jac2 		jac_ftheta_C1C2H4 				A1_C1C2H5 	A2_C1C2H5 	A1_C1C2H5 	D2_C1C2H5 	D1_C1C2H5 	jac1 	jac2 		jac_ftheta_C1C2H5 				A1_C1C2Br6 	A2_C1C2Br6 	A1_C1C2Br6 	D2_C1C2Br6 	D1_C1C2Br6 	jac1 	jac2 		jac_ftheta_C1C2Br6 				A1_C1H3H4 	A2_C1H3H4 	A1_C1H3H4 	D2_C1H3H4 	D1_C1H3H4 	jac1 	jac2 		jac_ftheta_C1H3H4 				A1_C1H3H5 	A2_C1H3H5 	A1_C1H3H5 	D2_C1H3H5 	D1_C1H3H5 	jac1 	jac2 		jac_ftheta_C1H3H5 				A1_C1H4H5 	A2_C1H4H5 	A1_C1H4H5 	D2_C1H4H5 	D1_C1H4H5 	jac1 	jac2 		jac_ftheta_C1H4H5 						A1_C2H3H4 	A2_C2H3H4 	A1_C2H3H4 	D2_C2H3H4 	D1_C2H3H4 	jac1 	jac2 		jac_ftheta_C2H3H4 				A1_C2H3H5 	A2_C2H3H5 	A1_C2H3H5 	D2_C2H3H5 	D1_C2H3H5 	jac1 	jac2 		jac_ftheta_C2H3H5 				A1_C2H4H5 	A2_C2H4H5 	A1_C2H4H5 	D2_C2H4H5 	D1_C2H4H5 	jac1 	jac2 		jac_ftheta_C2H4H5 					A1_C1H3Br6 	A2_C1H3Br6 	A1_C1H3Br6 	D2_C1H3Br6 	D1_C1H3Br6 	jac1 	jac2 		jac_ftheta_C1H3Br6 				A1_C1H4Br6 	A2_C1H4Br6 	A1_C1H4Br6 	D2_C1H4Br6 	D1_C1H4Br6 	jac1 	jac2 		jac_ftheta_C1H4Br6 					A1_C1H5Br6 	A2_C1H5Br6 	A1_C1H5Br6 	D2_C1H5Br6 	D1_C1H5Br6 	jac1 	jac2 		jac_ftheta_C1H5Br6 						A1_C2H3Br6 	A2_C2H3Br6 	A1_C2H3Br6 	D2_C2H3Br6 	D1_C2H3Br6 	jac1 	jac2 		jac_ftheta_C2H3Br6 				A1_C2H4Br6 	A2_C2H4Br6 	A1_C2H4Br6 	D2_C2H4Br6 	D1_C2H4Br6 	jac1 	jac2 		jac_ftheta_C2H4Br6 					A1_C2H5Br6 	A2_C2H5Br6 	A1_C2H5Br6 	D2_C2H5Br6 	D1_C2H5Br6 	jac1 	jac2 		jac_ftheta_C2H5Br6 					A1_H3H4H5 	A2_H3H4H5 	A1_H3H4H5 	D2_H3H4H5 	D1_H3H4H5 	jac1 	jac2 		jac_ftheta_H3H4H5 					A1_H3H4Br6 	A2_H3H4Br6 	A1_H3H4Br6 	D2_H3H4Br6 	D1_H3H4Br6 	jac1 	jac2 		jac_ftheta_H3H4Br6 					A1_H3H5Br6 	A2_H3H5Br6 	A1_H3H5Br6 	D2_H3H5Br6 	D1_H3H5Br6 	jac1 	jac2 		jac_ftheta_H3H5Br6 				A1_H4H5Br6 	A2_H4H5Br6 	A1_H4H5Br6 	D2_H4H5Br6 	D1_H4H5Br6 	jac1 	jac2 		jac_ftheta_H4H5Br6 


%--------------------------------------------------------------
%--------------------------------------------------------------
%--------------------------------------------------------------




%%%% START jac fPhi

W1HHHC = net_fphi_HHHC.IW{1,1};
W2HHHC = net_fphi_HHHC.LW{2,1};
B1HHHC = net_fphi_HHHC.b{1};
B2HHHC = net_fphi_HHHC.b{2};

W1HHHBr = net_fphi_HHHBr.IW{1,1};
W2HHHBr = net_fphi_HHHBr.LW{2,1};
B1HHHBr = net_fphi_HHHBr.b{1};
B2HHHBr = net_fphi_HHHBr.b{2};

W1HHCC = net_fphi_HHCC.IW{1,1};
W2HHCC = net_fphi_HHCC.LW{2,1};
B1HHCC = net_fphi_HHCC.b{1};
B2HHCC = net_fphi_HHCC.b{2};

W1HHCBr = net_fphi_HHCBr.IW{1,1};
W2HHCBr = net_fphi_HHCBr.LW{2,1};
B1HHCBr = net_fphi_HHCBr.b{1};
B2HHCBr = net_fphi_HHCBr.b{2};

W1HCCBr = net_fphi_HCCBr.IW{1,1};
W2HCCBr = net_fphi_HCCBr.LW{2,1};
B1HCCBr = net_fphi_HCCBr.b{1};
B2HCCBr = net_fphi_HCCBr.b{2};




S1 = 32;
S2 = 1;





A1_C1H3H4H5 = logsig(W1HHHC*[r_C1H3; r_C1H4 ; r_C1H5; r_H3H4; r_H3H5; r_H4H5] + (repmat(B1HHHC*ones(1,(1)*1),1,Q ) ));
A2_C1H3H4H5 = W2HHHC*A1_C1H3H4H5 + B2HHHC*ones(1,(1)*1);
A1_C1H3H4H5 = kron(A1_C1H3H4H5,ones(1,S2));
D2_C1H3H4H5 = nnmdlin(A2_C1H3H4H5);
D1_C1H3H4H5 = nnmdlog(A1_C1H3H4H5,D2_C1H3H4H5,W2HHHC);
jac1 = nnlmarq(kron([r_C1H3; r_C1H4 ; r_C1H5; r_H3H4; r_H3H5; r_H4H5],ones(1,S2)),D1_C1H3H4H5);
jac2 = nnlmarq(A1_C1H3H4H5,D2_C1H3H4H5);

jac_fphi_C1H3H4H5 = [jac1,D1_C1H3H4H5',jac2,D2_C1H3H4H5']; 

clear A1_C1H3H4H5 A2_C1H3H4H5 A1_C1H3H4H5 D2_C1H3H4H5 D1_C1H3H4H5 jac1 jac2 

% sim(net_fphi_HHHC,[r_C1H3; r_C1H4 ; r_C1H5; r_H3H4; r_H3H5; r_H4H5]) 





A1_C2H3H4H5 = logsig(W1HHHC*[r_C2H3; r_C2H4 ; r_C2H5; r_H3H4; r_H3H5; r_H4H5] + (repmat(B1HHHC*ones(1,(1)*1),1,Q ) ));
A2_C2H3H4H5 = W2HHHC*A1_C2H3H4H5 + B2HHHC*ones(1,(1)*1);
A1_C2H3H4H5 = kron(A1_C2H3H4H5,ones(1,S2));
D2_C2H3H4H5 = nnmdlin(A2_C2H3H4H5);
D1_C2H3H4H5 = nnmdlog(A1_C2H3H4H5,D2_C2H3H4H5,W2HHHC);
jac1 = nnlmarq(kron([r_C2H3; r_C2H4 ; r_C2H5; r_H3H4; r_H3H5; r_H4H5],ones(1,S2)),D1_C2H3H4H5);
jac2 = nnlmarq(A1_C2H3H4H5,D2_C2H3H4H5);

jac_fphi_C2H3H4H5 = [jac1,D1_C2H3H4H5',jac2,D2_C2H3H4H5']; 

clear A1_C2H3H4H5 A2_C2H3H4H5 A1_C2H3H4H5 D2_C2H3H4H5 D1_C2H3H4H5 jac1 jac2 

% sim(net_fphi_HHHC,[r_C2H3; r_C2H4 ; r_C2H5; r_H3H4; r_H3H5; r_H4H5])




A1_H3H4H5Br6 = logsig(W1HHHBr*[r_H3H4; r_H3H5 ; r_H3Br6; r_H4H5; r_H4Br6; r_H5Br6] + (repmat(B1HHHBr*ones(1,(1)*1),1,Q ) ));
A2_H3H4H5Br6 = W2HHHBr*A1_H3H4H5Br6 + B2HHHBr*ones(1,(1)*1);
A1_H3H4H5Br6 = kron(A1_H3H4H5Br6,ones(1,S2));
D2_H3H4H5Br6 = nnmdlin(A2_H3H4H5Br6);
D1_H3H4H5Br6 = nnmdlog(A1_H3H4H5Br6,D2_H3H4H5Br6,W2HHHBr);
jac1 = nnlmarq(kron([r_H3H4; r_H3H5 ; r_H3Br6; r_H4H5; r_H4Br6; r_H5Br6],ones(1,S2)),D1_H3H4H5Br6);
jac2 = nnlmarq(A1_H3H4H5Br6,D2_H3H4H5Br6);

jac_fphi_H3H4H5Br6 = [jac1,D1_H3H4H5Br6',jac2,D2_H3H4H5Br6']; 

clear A1_H3H4H5Br6 A2_H3H4H5Br6 A1_H3H4H5Br6 D2_H3H4H5Br6 D1_H3H4H5Br6 jac1 jac2    

% sim(net_fphi_HHHBr,[r_H3H4; r_H3H5 ; r_H3Br6; r_H4H5; r_H4Br6; r_H5Br6])




A1_C1C2H3H4 = logsig(W1HHCC*[r_C1C2; r_C1H3 ; r_C1H4; r_C2H3; r_C2H4; r_H3H4] + (repmat(B1HHCC*ones(1,(1)*1),1,Q ) ));
A2_C1C2H3H4 = W2HHCC*A1_C1C2H3H4 + B2HHCC*ones(1,(1)*1);
A1_C1C2H3H4 = kron(A1_C1C2H3H4,ones(1,S2));
D2_C1C2H3H4 = nnmdlin(A2_C1C2H3H4);
D1_C1C2H3H4 = nnmdlog(A1_C1C2H3H4,D2_C1C2H3H4,W2HHCC);
jac1 = nnlmarq(kron([r_C1C2; r_C1H3 ; r_C1H4; r_C2H3; r_C2H4; r_H3H4],ones(1,S2)),D1_C1C2H3H4);
jac2 = nnlmarq(A1_C1C2H3H4,D2_C1C2H3H4);

jac_fphi_C1C2H3H4 = [jac1,D1_C1C2H3H4',jac2,D2_C1C2H3H4']; 

clear A1_C1C2H3H4 A2_C1C2H3H4 A1_C1C2H3H4 D2_C1C2H3H4 D1_C1C2H3H4 jac1 jac2 

% sim(net_fphi_HHCC,[r_C1C2; r_C1H3 ; r_C1H4; r_C2H3; r_C2H4; r_H3H4])




A1_C1C2H3H5 = logsig(W1HHCC*[r_C1C2; r_C1H3 ; r_C1H5; r_C2H3; r_C2H5; r_H3H5] + (repmat(B1HHCC*ones(1,(1)*1),1,Q ) ));
A2_C1C2H3H5 = W2HHCC*A1_C1C2H3H5 + B2HHCC*ones(1,(1)*1);
A1_C1C2H3H5 = kron(A1_C1C2H3H5,ones(1,S2));
D2_C1C2H3H5 = nnmdlin(A2_C1C2H3H5);
D1_C1C2H3H5 = nnmdlog(A1_C1C2H3H5,D2_C1C2H3H5,W2HHCC);
jac1 = nnlmarq(kron([r_C1C2; r_C1H3 ; r_C1H5; r_C2H3; r_C2H5; r_H3H5],ones(1,S2)),D1_C1C2H3H5);
jac2 = nnlmarq(A1_C1C2H3H5,D2_C1C2H3H5);

jac_fphi_C1C2H3H5 = [jac1,D1_C1C2H3H5',jac2,D2_C1C2H3H5']; 

clear A1_C1C2H3H5 A2_C1C2H3H5 A1_C1C2H3H5 D2_C1C2H3H5 D1_C1C2H3H5 jac1 jac2 

% sim(net_fphi_HHCC,[r_C1C2; r_C1H3 ; r_C1H5; r_C2H3; r_C2H5; r_H3H5])




A1_C1C2H4H5 = logsig(W1HHCC*[r_C1C2; r_C1H4 ; r_C1H5; r_C2H4; r_C2H5; r_H4H5] + (repmat(B1HHCC*ones(1,(1)*1),1,Q ) ));
A2_C1C2H4H5 = W2HHCC*A1_C1C2H4H5 + B2HHCC*ones(1,(1)*1);
A1_C1C2H4H5 = kron(A1_C1C2H4H5,ones(1,S2));
D2_C1C2H4H5 = nnmdlin(A2_C1C2H4H5);
D1_C1C2H4H5 = nnmdlog(A1_C1C2H4H5,D2_C1C2H4H5,W2HHCC);
jac1 = nnlmarq(kron([r_C1C2; r_C1H4 ; r_C1H5; r_C2H4; r_C2H5; r_H4H5],ones(1,S2)),D1_C1C2H4H5);
jac2 = nnlmarq(A1_C1C2H4H5,D2_C1C2H4H5);

jac_fphi_C1C2H4H5 = [jac1,D1_C1C2H4H5',jac2,D2_C1C2H4H5']; 

clear A1_C1C2H4H5 A2_C1C2H4H5 A1_C1C2H4H5 D2_C1C2H4H5 D1_C1C2H4H5 jac1 jac2 

% sim(net_fphi_HHCC,[r_C1C2; r_C1H4 ; r_C1H5; r_C2H4; r_C2H5; r_H4H5])





A1_C1H3H4Br6 = logsig(W1HHCBr*[r_C1H3; r_C1H4 ; r_C1Br6; r_H3H4; r_H3Br6; r_H4Br6] + (repmat(B1HHCBr*ones(1,(1)*1),1,Q ) ));
A2_C1H3H4Br6 = W2HHCBr*A1_C1H3H4Br6 + B2HHCBr*ones(1,(1)*1);
A1_C1H3H4Br6 = kron(A1_C1H3H4Br6,ones(1,S2));
D2_C1H3H4Br6 = nnmdlin(A2_C1H3H4Br6);
D1_C1H3H4Br6 = nnmdlog(A1_C1H3H4Br6,D2_C1H3H4Br6,W2HHCBr);
jac1 = nnlmarq(kron([r_C1H3; r_C1H4 ; r_C1Br6; r_H3H4; r_H3Br6; r_H4Br6],ones(1,S2)),D1_C1H3H4Br6);
jac2 = nnlmarq(A1_C1H3H4Br6,D2_C1H3H4Br6);

jac_fphi_C1H3H4Br6 = [jac1,D1_C1H3H4Br6',jac2,D2_C1H3H4Br6']; 

clear A1_C1H3H4Br6 A2_C1H3H4Br6 A1_C1H3H4Br6 D2_C1H3H4Br6 D1_C1H3H4Br6 jac1 jac2 

% sim(net_fphi_HHCBr,[r_C1H3; r_C1H4 ; r_C1Br6; r_H3H4; r_H3Br6; r_H4Br6])




A1_C1H3H5Br6 = logsig(W1HHCBr*[r_C1H3; r_C1H5 ; r_C1Br6; r_H3H5; r_H3Br6; r_H5Br6] + (repmat(B1HHCBr*ones(1,(1)*1),1,Q ) ));
A2_C1H3H5Br6 = W2HHCBr*A1_C1H3H5Br6 + B2HHCBr*ones(1,(1)*1);
A1_C1H3H5Br6 = kron(A1_C1H3H5Br6,ones(1,S2));
D2_C1H3H5Br6 = nnmdlin(A2_C1H3H5Br6);
D1_C1H3H5Br6 = nnmdlog(A1_C1H3H5Br6,D2_C1H3H5Br6,W2HHCBr);
jac1 = nnlmarq(kron([r_C1H3; r_C1H5 ; r_C1Br6; r_H3H5; r_H3Br6; r_H5Br6],ones(1,S2)),D1_C1H3H5Br6);
jac2 = nnlmarq(A1_C1H3H5Br6,D2_C1H3H5Br6);

jac_fphi_C1H3H5Br6 = [jac1,D1_C1H3H5Br6',jac2,D2_C1H3H5Br6']; 

clear A1_C1H3H5Br6 A2_C1H3H5Br6 A1_C1H3H5Br6 D2_C1H3H5Br6 D1_C1H3H5Br6 jac1 jac2 

% sim(net_fphi_HHCBr,[r_C1H3; r_C1H5 ; r_C1Br6; r_H3H5; r_H3Br6; r_H5Br6])




A1_C1H4H5Br6 = logsig(W1HHCBr*[r_C1H4; r_C1H5 ; r_C1Br6; r_H4H5; r_H4Br6; r_H5Br6] + (repmat(B1HHCBr*ones(1,(1)*1),1,Q ) ));
A2_C1H4H5Br6 = W2HHCBr*A1_C1H4H5Br6 + B2HHCBr*ones(1,(1)*1);
A1_C1H4H5Br6 = kron(A1_C1H4H5Br6,ones(1,S2));
D2_C1H4H5Br6 = nnmdlin(A2_C1H4H5Br6);
D1_C1H4H5Br6 = nnmdlog(A1_C1H4H5Br6,D2_C1H4H5Br6,W2HHCBr);
jac1 = nnlmarq(kron([r_C1H4; r_C1H5 ; r_C1Br6; r_H4H5; r_H4Br6; r_H5Br6],ones(1,S2)),D1_C1H4H5Br6);
jac2 = nnlmarq(A1_C1H4H5Br6,D2_C1H4H5Br6);

jac_fphi_C1H4H5Br6 = [jac1,D1_C1H4H5Br6',jac2,D2_C1H4H5Br6']; 

clear A1_C1H4H5Br6 A2_C1H4H5Br6 A1_C1H4H5Br6 D2_C1H4H5Br6 D1_C1H4H5Br6 jac1 jac2 

% sim(net_fphi_HHCBr,[r_C1H4; r_C1H5 ; r_C1Br6; r_H4H5; r_H4Br6; r_H5Br6])






A1_C2H3H4Br6 = logsig(W1HHCBr*[r_C2H3; r_C2H4 ; r_C2Br6; r_H3H4; r_H3Br6; r_H4Br6] + (repmat(B1HHCBr*ones(1,(1)*1),1,Q ) ));
A2_C2H3H4Br6 = W2HHCBr*A1_C2H3H4Br6 + B2HHCBr*ones(1,(1)*1);
A1_C2H3H4Br6 = kron(A1_C2H3H4Br6,ones(1,S2));
D2_C2H3H4Br6 = nnmdlin(A2_C2H3H4Br6);
D1_C2H3H4Br6 = nnmdlog(A1_C2H3H4Br6,D2_C2H3H4Br6,W2HHCBr);
jac1 = nnlmarq(kron([r_C2H3; r_C2H4 ; r_C2Br6; r_H3H4; r_H3Br6; r_H4Br6],ones(1,S2)),D1_C2H3H4Br6);
jac2 = nnlmarq(A1_C2H3H4Br6,D2_C2H3H4Br6);

jac_fphi_C2H3H4Br6 = [jac1,D1_C2H3H4Br6',jac2,D2_C2H3H4Br6']; 

clear A1_C2H3H4Br6 A2_C2H3H4Br6 A1_C2H3H4Br6 D2_C2H3H4Br6 D1_C2H3H4Br6 jac1 jac2 

% sim(net_fphi_HHCBr,[r_C2H3; r_C2H4 ; r_C2Br6; r_H3H4; r_H3Br6; r_H4Br6])




A1_C2H3H5Br6 = logsig(W1HHCBr*[r_C2H3; r_C2H5 ; r_C2Br6; r_H3H5; r_H3Br6; r_H5Br6] + (repmat(B1HHCBr*ones(1,(1)*1),1,Q ) ));
A2_C2H3H5Br6 = W2HHCBr*A1_C2H3H5Br6 + B2HHCBr*ones(1,(1)*1);
A1_C2H3H5Br6 = kron(A1_C2H3H5Br6,ones(1,S2));
D2_C2H3H5Br6 = nnmdlin(A2_C2H3H5Br6);
D1_C2H3H5Br6 = nnmdlog(A1_C2H3H5Br6,D2_C2H3H5Br6,W2HHCBr);
jac1 = nnlmarq(kron([r_C2H3; r_C2H5 ; r_C2Br6; r_H3H5; r_H3Br6; r_H5Br6],ones(1,S2)),D1_C2H3H5Br6);
jac2 = nnlmarq(A1_C2H3H5Br6,D2_C2H3H5Br6);

jac_fphi_C2H3H5Br6 = [jac1,D1_C2H3H5Br6',jac2,D2_C2H3H5Br6']; 

clear A1_C2H3H5Br6 A2_C2H3H5Br6 A1_C2H3H5Br6 D2_C2H3H5Br6 D1_C2H3H5Br6 jac1 jac2 

% sim(net_fphi_HHCBr,[r_C2H3; r_C2H5 ; r_C2Br6; r_H3H5; r_H3Br6; r_H5Br6])




A1_C2H4H5Br6 = logsig(W1HHCBr*[r_C2H4; r_C2H5 ; r_C2Br6; r_H4H5; r_H4Br6; r_H5Br6] + (repmat(B1HHCBr*ones(1,(1)*1),1,Q ) ));
A2_C2H4H5Br6 = W2HHCBr*A1_C2H4H5Br6 + B2HHCBr*ones(1,(1)*1);
A1_C2H4H5Br6 = kron(A1_C2H4H5Br6,ones(1,S2));
D2_C2H4H5Br6 = nnmdlin(A2_C2H4H5Br6);
D1_C2H4H5Br6 = nnmdlog(A1_C2H4H5Br6,D2_C2H4H5Br6,W2HHCBr);
jac1 = nnlmarq(kron([r_C2H4; r_C2H5 ; r_C2Br6; r_H4H5; r_H4Br6; r_H5Br6],ones(1,S2)),D1_C2H4H5Br6);
jac2 = nnlmarq(A1_C2H4H5Br6,D2_C2H4H5Br6);

jac_fphi_C2H4H5Br6 = [jac1,D1_C2H4H5Br6',jac2,D2_C2H4H5Br6']; 

clear A1_C2H4H5Br6 A2_C2H4H5Br6 A1_C2H4H5Br6 D2_C2H4H5Br6 D1_C2H4H5Br6 jac1 jac2 

% sim(net_fphi_HHCBr,[r_C2H4; r_C2H5 ; r_C2Br6; r_H4H5; r_H4Br6; r_H5Br6])





A1_C1C2H3Br6 = logsig(W1HCCBr*[r_C1C2; r_C1H3 ; r_C1Br6; r_C2H3; r_C2Br6; r_H3Br6] + (repmat(B1HCCBr*ones(1,(1)*1),1,Q ) ));
A2_C1C2H3Br6 = W2HCCBr*A1_C1C2H3Br6 + B2HCCBr*ones(1,(1)*1);
A1_C1C2H3Br6 = kron(A1_C1C2H3Br6,ones(1,S2));
D2_C1C2H3Br6 = nnmdlin(A2_C1C2H3Br6);
D1_C1C2H3Br6 = nnmdlog(A1_C1C2H3Br6,D2_C1C2H3Br6,W2HCCBr);
jac1 = nnlmarq(kron([r_C1C2; r_C1H3 ; r_C1Br6; r_C2H3; r_C2Br6; r_H3Br6],ones(1,S2)),D1_C1C2H3Br6);
jac2 = nnlmarq(A1_C1C2H3Br6,D2_C1C2H3Br6);

jac_fphi_C1C2H3Br6 = [jac1,D1_C1C2H3Br6',jac2,D2_C1C2H3Br6']; 

clear A1_C1C2H3Br6 A2_C1C2H3Br6 A1_C1C2H3Br6 D2_C1C2H3Br6 D1_C1C2H3Br6 jac1 jac2 

% sim(net_fphi_HCCBr,[r_C1C2; r_C1H3 ; r_C1Br6; r_C2H3; r_C2Br6; r_H3Br6])





A1_C1C2H4Br6 = logsig(W1HCCBr*[r_C1C2; r_C1H4 ; r_C1Br6; r_C2H4; r_C2Br6; r_H4Br6] + (repmat(B1HCCBr*ones(1,(1)*1),1,Q ) ));
A2_C1C2H4Br6 = W2HCCBr*A1_C1C2H4Br6 + B2HCCBr*ones(1,(1)*1);
A1_C1C2H4Br6 = kron(A1_C1C2H4Br6,ones(1,S2));
D2_C1C2H4Br6 = nnmdlin(A2_C1C2H4Br6);
D1_C1C2H4Br6 = nnmdlog(A1_C1C2H4Br6,D2_C1C2H4Br6,W2HCCBr);
jac1 = nnlmarq(kron([r_C1C2; r_C1H4 ; r_C1Br6; r_C2H4; r_C2Br6; r_H4Br6],ones(1,S2)),D1_C1C2H4Br6);
jac2 = nnlmarq(A1_C1C2H4Br6,D2_C1C2H4Br6);

jac_fphi_C1C2H4Br6 = [jac1,D1_C1C2H4Br6',jac2,D2_C1C2H4Br6']; 

clear A1_C1C2H4Br6 A2_C1C2H4Br6 A1_C1C2H4Br6 D2_C1C2H4Br6 D1_C1C2H4Br6 jac1 jac2 

% sim(net_fphi_HCCBr,[r_C1C2; r_C1H4 ; r_C1Br6; r_C2H4; r_C2Br6; r_H4Br6])




A1_C1C2H5Br6 = logsig(W1HCCBr*[r_C1C2; r_C1H5 ; r_C1Br6; r_C2H5; r_C2Br6; r_H5Br6] + (repmat(B1HCCBr*ones(1,(1)*1),1,Q ) ));
A2_C1C2H5Br6 = W2HCCBr*A1_C1C2H5Br6 + B2HCCBr*ones(1,(1)*1);
A1_C1C2H5Br6 = kron(A1_C1C2H5Br6,ones(1,S2));
D2_C1C2H5Br6 = nnmdlin(A2_C1C2H5Br6);
D1_C1C2H5Br6 = nnmdlog(A1_C1C2H5Br6,D2_C1C2H5Br6,W2HCCBr);
jac1 = nnlmarq(kron([r_C1C2; r_C1H5 ; r_C1Br6; r_C2H5; r_C2Br6; r_H5Br6],ones(1,S2)),D1_C1C2H5Br6);
jac2 = nnlmarq(A1_C1C2H5Br6,D2_C1C2H5Br6);

jac_fphi_C1C2H5Br6 = [jac1,D1_C1C2H5Br6',jac2,D2_C1C2H5Br6']; 

clear A1_C1C2H5Br6 A2_C1C2H5Br6 A1_C1C2H5Br6 D2_C1C2H5Br6 D1_C1C2H5Br6 jac1 jac2 

% sim(net_fphi_HCCBr,[r_C1C2; r_C1H5 ; r_C1Br6; r_C2H5; r_C2Br6; r_H5Br6]);




jac_fphi_HHHC= [jac_fphi_C1H3H4H5 + jac_fphi_C2H3H4H5];
jac_fphi_HHHBr= [jac_fphi_H3H4H5Br6];
jac_fphi_HHCC= [jac_fphi_C1C2H3H4 + jac_fphi_C1C2H3H5 + jac_fphi_C1C2H4H5];
jac_fphi_HHCBr= [jac_fphi_C1H3H4Br6 + jac_fphi_C1H3H5Br6 + jac_fphi_C1H4H5Br6 + jac_fphi_C2H3H4Br6 + jac_fphi_C2H3H5Br6 + jac_fphi_C2H4H5Br6];
jac_fphi_HCCBr= [jac_fphi_C1C2H3Br6 + jac_fphi_C1C2H4Br6 + jac_fphi_C1C2H5Br6];



clear jac_fphi_C1H3H4H5 jac_fphi_C2H3H4H5 jac_fphi_H3H4H5Br6 jac_fphi_C1C2H3H4 jac_fphi_C1C2H3H5 jac_fphi_C1C2H4H5;
clear jac_fphi_C1H3H4Br6 jac_fphi_C1H3H5Br6 jac_fphi_C1H4H5Br6 jac_fphi_C2H3H4Br6 jac_fphi_C2H3H5Br6 jac_fphi_C2H4H5Br6;
clear jac_fphi_C1C2H3Br6 jac_fphi_C1C2H4Br6 jac_fphi_C1C2H5Br6

J_fphi_all= [jac_fphi_HHHC jac_fphi_HHHBr jac_fphi_HHCC jac_fphi_HHCBr jac_fphi_HCCBr];

clear jac_fphi_HHHC jac_fphi_HHHBr jac_fphi_HHCC jac_fphi_HHCBr jac_fphi_HCCBr
clear W1HHHC 	W2HHHC 	B1HHHC 	B2HHHC 		W1HHHBr 	W2HHHBr 	B1HHHBr 	B2HHHBr 		W1HHCC 	W2HHCC 	B1HHCC 	B2HHCC 		W1HHCBr 	W2HHCBr 	B1HHCBr 	B2HHCBr 		W1HCCBr 	W2HCCBr 	B1HCCBr 	B2HCCBr 


%%%% END jac fPhi


clear clust_size net_fc net_fphi_HCCBr net_fphi_HHCBr net_fphi_HHCC net_fphi_HHHBr net_fphi_HHHC net_fr_CBr net_fr_CC net_fr_HBr net_fr_HC net_fr_HH
clear net_ftheta_CCBr net_ftheta_HCBr net_ftheta_HCC net_ftheta_HHBr net_ftheta_HHC net_ftheta_HHH

clear r_C1Br6 r_C1C2 r_C1H3 r_C1H4 r_C1H5 r_C2Br6 r_C2H3 r_C2H4 r_C2H5 r_H3Br6 r_H3H4 r_H3H5 r_H4Br6 r_H4H5 r_H5Br6;
clear clust_size;




J_Milind =horzcat(J_fr_all,J_ftheta_all,J_fphi_all);

clear J_fr_all J_ftheta_all J_fphi_all

JE_Milind= J_Milind'*Ex;
JtJ_Milind= J_Milind'*J_Milind;
normJE_Milind= sqrt(JE_Milind'*JE_Milind);
%%%----------------------------------------------------

i=1; %dummy line for breakpoint


% jac_ftheta_H1H2Br = [jac1,D1_H1H2Br',jac2,D2_H1H2Br']; 
% 
% J_ftheta_H1H2Br= jac_ftheta_H1H2Br;
% 
% 
% J_Milind =horzcat(J_fr_HH_HBr,J_ftheta_H1H2Br);
% JE_Milind= J_Milind'*Ex';
% JtJ_Milind= J_Milind'*J_Milind;
% normJE_Milind= sqrt(JE_Milind'*JE_Milind);

%%%%%$$$$$$$$$$@@@@@






% % %%%%%%%
% S1 = 40;
% S2 = 1;
% for iQ = 1:Q
% 
%     ciQ = 0;
%     ciQHH = 0; ciQHBr = 0; 
% 
%     for i=1:1:clust_size(iQ)
%         for j=1:1:clust_size(iQ)
%             if j>i%j~= i %
%                 ciQ = ciQ+1;
%                 identify2=sprintf('%d %d', type(i),type(j));
%                 switch identify2
%                    case {'1 1'}
%                       ciQHH = ciQHH+1;
%                       r_fr_HH(ciQHH) = rs(i,j,iQ);                      
%                    case {'1 35', '35 1'}
%                       ciQHBr = ciQHBr+1;
%                       r_fr_HBr(ciQHBr) = rs(i,j,iQ);
%                    
%                 end
%                 r_fr(ciQ) = rs(i,j,iQ);
% 
%             end
%         end
%     end
% 
% 
% 
% 
%     A1HH = logsig(W1HH*r_fr_HH + B1HH*ones(1,(1)*1));
%     A2HH = W2HH*A1HH + B2HH*ones(1,(1)*1);
%     % FIND JACOBIAN
%     A1HH = kron(A1HH,ones(1,S2));
%     D2HH = nnmdlin(A2HH);
% 
%     D1HH = nnmdlog(A1HH,D2HH,W2HH);
%     jac1 = nnlmarq(kron(r_fr_HH,ones(1,S2)),D1HH);
%     jac2 = nnlmarq(A1HH,D2HH);
% 
%     jac_fr_HH = [jac1,D1HH',jac2,D2HH']; % Make jacbian for current iQ configuration
%     
%     for j= 1: 3*S1+1
%         jac(1,j)= 0;
%         for i = 1:(1)
%             jac(1,j) = jac(1,j)+jac_fr_HH(i,j);  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
%     
%     
%     A1HBr = logsig(W1HBr*r_fr_HBr + B1HBr*ones(1,(2)*1));
%     A2HBr = W2HBr*A1HBr + B2HBr*ones(1,(2)*1);
%     % FIND JACOBIAN
%     A1HBr = kron(A1HBr,ones(1,S2));
%     D2HBr = nnmdlin(A2HBr);
% 
%     D1HBr = nnmdlog(A1HBr,D2HBr,W2HBr);
%     jac1 = nnlmarq(kron(r_fr_HBr,ones(1,S2)),D1HBr);
%     jac2 = nnlmarq(A1HBr,D2HBr);
% 
%     jac_fr_HBr = [jac1,D1HBr',jac2,D2HBr']; % Make jacbian for current iQ configuration
% 
%     for j= (3*S1+1)+1: 2*(3*S1+1)
%         jac(1,j)= 0;
%         for i = 1:(2)
%             jac(1,j) = jac(1,j)+jac_fr_HBr(i,j-(3*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
%     
%     for j= 1:2*(3*S1+1)
%         J_fr(iQ,j)= jac(1,j);
%     end
% 
% end
% 
% 
% % %%%%%%%
% 
% 
% % J_theta will have dimensions of Q x ( 5N + 1)
% % First layer W (N x 3), b (N x 1), Second Layer W (1 x N ), b (1x1)
% 
% 
% 
% W1HHBr = net_ftheta_HHBr.IW{1,1};
% W2HHBr = net_ftheta_HHBr.LW{2,1};
% B1HHBr = net_ftheta_HHBr.b{1};
% B2HHBr = net_ftheta_HHBr.b{2};
% 
% S1 = 120;
% S2 = 1;
% 
% 
% 
% 
% for iQ = 1:Q
%     ciQ = 0;
%     ciQHHBr=0;
%     for i=1:1:clust_size(iQ)
%         for j=1:1:clust_size(iQ)
%             for k=1:1:clust_size(iQ)
%                 if (j>i && k>i && k>j)
%                     ciQ = ciQ+1;
%                     identify3=sprintf('%d %d %d', type(i),type(j),type(k));
%                     switch identify3
%                                               
%                       case {'1 1 35','1 35 1','35 1 1'}
%                            ciQHHBr=ciQHHBr+1;
%                            r_ftheta_HHBr(1,ciQHHBr)=rs(i,j,iQ);
%                            r_ftheta_HHBr(2,ciQHHBr)=rs(i,k,iQ);
%                            r_ftheta_HHBr(3,ciQHHBr)=rs(j,k,iQ);                        
%                       
%                     end
% %                     r_ftheta(1,ciQ) = rs(i,j,iQ);
% %                     r_ftheta(2,ciQ) = rs(i,k,iQ);
% %                     r_ftheta(3,ciQ) = rs(j,k,iQ);
% 
%                 end
%             end
%         end
%     end
% 
%     
%     A1HHBr = logsig(W1HHBr*r_ftheta_HHBr + B1HHBr*ones(1,(1)*1));
%     A2HHBr = W2HHBr*A1HHBr+B2HHBr*ones(1,(1)*1);
% 
%     % FIND JACOBIAN
%     A1HHBr = kron(A1HHBr,ones(1,S2));
%     D2HHBr = nnmdlin(A2HHBr);
% 
%     D1HHBr = nnmdlog(A1HHBr,D2HHBr,W2HHBr);
%     jac1 = nnlmarq(kron(r_ftheta_HHBr,ones(1,S2)),D1HHBr);
%     jac2 = nnlmarq(A1HHBr,D2HHBr);
% 
%     jac_ftheta_HHBr = [jac1,D1HHBr',jac2,D2HHBr'];
%     %%%%
%     for j= 0*(5*S1+1)+1: 1*(5*S1+1)
%         jac(1,j)= 0;
%         for i = 1:(1)
%             jac(1,j) = jac(1,j)+jac_ftheta_HHBr(i,j-0*(5*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
%     
%     for j= 1: 1*(5*S1+1)
%         J_ftheta(iQ,j)= jac(1,j);
%     end
%     
% end
% 
%  
% 
% 
% 
% 
% %%%Make cobined Jacobian matrix
% % J = zeros(Q,5*(3*S1+1)+6*(5*S1+1)+5*(8*S1+1) );
% 
% % J(:,1:5*(3*S1+1))=J_fr;
% % J(:,5*(3*S1+1)+1:5*(3*S1+1)+6*(5*S1+1))=J_ftheta;
% % J(:,(3*S1+1)+(5*S1+1)+1:(3*S1+1)+(5*S1+1)+(8*S1+1))=J_fphi;
% J=horzcat(J_fr,J_ftheta);
% JE=J'*Ex';
% JtJ=J'*J;
% normJE=sqrt(JE'*JE);



% i=1; %dummy line for breakpoint





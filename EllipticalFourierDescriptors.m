(* ::Package:: *)

(* ::Section:: *)
(*Begin Package*)


BeginPackage["EllipticalDescriptors`"]


generateContour::usage = "generateContour[input] either takes in a binarized image or a list of positions and creates a contour";
ellipticalCoefficients::usage = "elliptical Coefficients[contour, order] takes in a contour and a harmonic order to compute elliptical
fourier coefficients";
plotEFD::usage = "plotEFD[contour_, coeff_, ptsToplot_:300, locus_:None] takes a contour and the elliptical fourier coefficients
to plot the shape from the coefficients atop the contour";
lobeContribution::usage = "lobeContribution[contour,ellipticalcoefficients] will take the contour and elliptical coefficients to give
the lobe contribution effects. The metric is based on the work of YE Sanchez-Corrales et al";
EllipticalFourierDescriptor::usage = "EllipticalFourierDescriptors[contour, options] takes in a binarized image or a list of positions
to compute the elliptical fourier coefficients. Different options for which the default values are set as follows:
\"normalize\" -> True, \"sizeInvariance\" -> True, \"order\" -> 10, \"ptsToplot\" ->300, \"locus\" -> {0.,0.}, \"lobeContribution\" ->
False";


(* ::Section:: *)
(*Main Function*)


Begin["`Private`"];


generateContour[input_]:=Module[{contourpositions,shortestpath},
contourpositions=Switch[Head@input,Image,
DeleteDuplicates@PixelValuePositions[MorphologicalPerimeter[input],1],_,
DeleteDuplicates@input];
shortestpath=Last@FindShortestTour[contourpositions];
Part[contourpositions,shortestpath]
];


paramsContour[contour_]:=Module[{dXY,dt,t,T,\[Phi]},
dXY=Differences[contour];
dt=(x\[Function]N@Sqrt[x])[Total/@(dXY^2)];
t=Prepend[Accumulate[dt],0];
T=Last@t;
\[Phi]=2 \[Pi] t/T;
{dXY,dt,t,T,\[Phi]}
];


ellipticalCoefficients[contour_,order_]:=Module[{computeCoefficients,coefficientsL,dt,t,dXY,\[Phi],T},
{dXY,dt,t,T,\[Phi]}=paramsContour[contour];
computeCoefficients[n_Integer]:=Block[{const=T/(2 (n^2) (\[Pi]^2)),
phiN=\[Phi]*n,dCosPhiN,dSinPhiN,an,bn,cn,dn},
dCosPhiN=Cos[Rest@phiN]-Cos[Most@phiN];
dSinPhiN=Sin[Rest@phiN]-Sin[Most@phiN];
an=const*Total@(dXY[[All,1]]/dt*dCosPhiN);
bn=const*Total@(dXY[[All,1]]/dt*dSinPhiN);
cn=const*Total@(dXY[[All,2]]/dt*dCosPhiN);
dn=const*Total@(dXY[[All,2]]/dt*dSinPhiN);
{an,bn,cn,dn}
];
computeCoefficients/@Range[order]
];


normalizeEFD[coeff_,sizeInvariance_:True]:=With[{\[Theta]=0.5 ArcTan[(coeff[[1,1]]^2-coeff[[1,2]]^2+coeff[[1,3]]^2-coeff[[1,4]]^2),
2 (coeff[[1,1]] coeff[[1,2]]+coeff[[1,3]] coeff[[1,4]])]},
Block[{coeffList,\[Psi],\[Psi]r,a,c,\[Theta]r},
\[Theta]r[i_Integer]:={{Cos[i \[Theta]],-Sin[i \[Theta]]},{Sin[i \[Theta]],Cos[i \[Theta]]}};
coeffList=Partition[#,2,2]&/@coeff;
coeffList=MapIndexed[(#1.\[Theta]r[First@#2])&,coeffList];
a=First@*Flatten@*First@coeffList;
c=(Flatten@First@coeffList)[[3]];
\[Psi]=ArcTan[a,c];
\[Psi]r={{Cos[\[Psi]],Sin[\[Psi]]},{-Sin[\[Psi]],Cos[\[Psi]]}};
coeffList=MapIndexed[Flatten[\[Psi]r.#]&,coeffList];
If[sizeInvariance,a=First@@coeffList;coeffList/=a,coeffList]
]
];


plotEFD[contour_,coeff_,ptsToplot_:300,locus_:None]:=First@Last@Reap@Block[{order=Length@coeff,nhalf,ti,xt,yt,\[Theta],dcVal,
dcCoefficients, xi,A0,\[Delta],C0,dXY,dt,t,T,\[Phi]},
dcCoefficients:=({dXY,dt,t,T,\[Phi]}=paramsContour[contour];
xi=Accumulate[Part[dXY,All,1]]-(dXY[[All,1]]/dt)*Rest@t;
A0=(1/T)*Total@(dXY[[All,1]]/(2 dt)*Differences[t^2]+xi*dt);
\[Delta]=Accumulate[Part[dXY,All,2]]-(dXY[[All,2]]/dt)*Rest@t;
C0=(1/T)*Total@(dXY[[All,2]]/(2 dt)*Differences[t^2]+\[Delta]*dt);
{First@First@contour+A0,Last@First@contour+C0});
nhalf=Ceiling[order/2];
ti=Array[#&,ptsToplot,{0.,1}];
dcVal=If[locus==None,dcCoefficients,locus];
xt=ConstantArray[dcVal[[1]],ptsToplot];
yt=ConstantArray[dcVal[[2]],ptsToplot];
Do[
xt+=coeff[[j,1]] Cos[2 (j) \[Pi] ti]+coeff[[j,2]] Sin[2 (j) \[Pi] ti];
yt+=coeff[[j,3]] Cos[2 (j) \[Pi] ti]+coeff[[j,4]] Sin[2 (j) \[Pi] ti];
Sow[
Graphics[{{Thick,XYZColor[0,0,1,0.6],Line@Thread[{xt,yt}]},{Thick,Dashed,Black,Line@contour},
Inset[Style[#,{XYZColor[1,0,0,0.8],Bold,FontSize->12}]&@Text["harmonic: "<>ToString@j]]}]
],{j,order}]
];


lobeContribution[contour:{{__Integer}..},modes_]:=Module[{starmat={},lcL={},lambdaminus,lambdaminusvec={},
lambdaplus,lambdaplusvec={},zetaplus,zetaplusvec={},zetaminus,zetaminusvec={},lczetaplus,lclambdaplus,lczetaminus,lclambdaminus,
lcaplus,lcbplus,lccplus,lcdplus,lcaminus,lcbminus,lccminus,lcdminus,maxerror=1 10^-13,tau,alphaprime,gammaprime,\[Rho],alphastar,
betastar,gammastar,deltastar,r,a,b,c,d,\[Phi],aprime,bprime,cprime,dprime,\[Theta],lambda1,lambda12,lambda21,lambda2,
lclambdaplusZeroM,lczetaplusZeroM,len=Length@modes},

lczetaplus=lclambdaplus=lczetaminus=lclambdaminus=ConstantArray[0,len];
lcaplus=lcbplus=lccplus=lcdplus=lcaminus=lcbminus=lccminus=lcdminus=ConstantArray[0,len];
tau=0.5 ArcTan[modes[[1,1]]^2+modes[[1,3]]^2-modes[[1,2]]^2-modes[[1,-1]]^2,
2.0*(modes[[1,1]] modes[[1,2]]+modes[[1,3]] modes[[1,-1]])];
alphaprime=modes[[1,1]] Cos[tau]+modes[[1,2]] Sin[tau];
gammaprime=modes[[1,3]] Cos[tau]+modes[[1,-1]] Sin[tau];
\[Rho]=ArcTan[alphaprime,gammaprime];
If[\[Rho]<0,tau+=N[\[Pi]]];

starmat=First@Last@Reap[Do[
alphastar=modes[[k,1]] Cos[k tau]+modes[[k,2]] Sin[k tau];
betastar=-modes[[k,1]] Sin[k tau]+modes[[k,2]] Cos[k tau];
gammastar=modes[[k,3]] Cos[k tau]+modes[[k,-1]] Sin[k tau];
deltastar=-modes[[k,3]] Sin[k tau]+modes[[k,-1]] Cos[k tau];
Sow[{alphastar,betastar,gammastar,deltastar}];
,{k,len}]];

r=starmat[[1,1]] starmat[[1,-1]]-starmat[[1,2]] starmat[[1,-2]];
If[r<0,starmat[[All,{2,-1}]]*=-1];
{a,b,c,d}=starmat\[Transpose];

Do[
\[Phi]=0.5 ArcTan[a[[i]]^2+c[[i]]^2-b[[i]]^2-d[[i]]^2,2.0 (a[[i]] b[[i]]+c[[i]] d[[i]])];
aprime=a[[i]] Cos[\[Phi]]+b[[i]] Sin[\[Phi]];
bprime=-a[[i]] Sin[\[Phi]]+b[[i]] Cos[\[Phi]];
cprime=c[[i]] Cos[\[Phi]]+d[[i]] Sin[\[Phi]];
dprime=-c[[i]] Sin[\[Phi]]+d[[i]] Cos[\[Phi]];
\[Theta]=ArcTan[aprime,cprime];
lambda1=aprime Cos[\[Theta]]+cprime Sin[\[Theta]];
lambda12=bprime Cos[\[Theta]]+dprime Sin[\[Theta]];
lambda21=-aprime Sin[\[Theta]]+cprime Cos[\[Theta]];
lambda2=-bprime Sin[\[Theta]]+dprime Cos[\[Theta]];
If[Or[Abs[lambda12]>maxerror,Abs[lambda21]>maxerror,lambda1<0,lambda1<Abs@lambda2],
Which[Abs[lambda12]>maxerror||Abs[lambda21]>maxerror,
Print["off-diagonal lambda matrix unequal to zero:\n"],lambda1<0,Print["Warning: lambda1 negative\n"],
lambda1<Abs[lambda2],Print["Warning: lambda 1 < |lambda2|:\n"];
Print["mode:\n","(lambda1,lambda12): ",{lambda1,lambda12},"\n","(lambda21,lambda2): ",{lambda1,lambda12},"\n"];];];
lambdaplus=(lambda1+lambda2)/2.0;
lambdaminus=(lambda1-lambda2)/2.0;
zetaplus=\[Theta]-\[Phi];
zetaminus=-\[Theta]-\[Phi];
AppendTo[zetaplusvec,zetaplus];
AppendTo[lambdaplusvec,lambdaplus];
AppendTo[lambdaminusvec,lambdaminus];
AppendTo[zetaminusvec,zetaminus],
{i,Length@starmat}];

lclambdaplusZeroM=lambdaplusvec[[2]];
lczetaplusZeroM=zetaplusvec[[2]];
lclambdaplus[[1]]=lambdaplusvec[[1]];
lczetaplus[[1]]=zetaplusvec[[1]];

Do[lclambdaplus[[j]]=lambdaplusvec[[j+1]];
lczetaplus[[j]]=zetaplusvec[[j+1]],{j,2,len-1}];

Do[lclambdaminus[[j]]=lambdaminusvec[[j-1]];
lczetaminus[[j]]=zetaminusvec[[j-1]],{j,2,len}];

If[True,Print["Ln quadruplets:\n ======================"];
Print["lc-EFA mode 0: ",{lclambdaplusZeroM,lczetaplusZeroM}];
Do[Print["lc-EFA mode ",i,"\n","(lambdaplus: ",lclambdaplus[[i]]," lambdaminus: ",lclambdaminus[[i]],")\n","(zetaplus: ",
lczetaplus[[i]]," zetaminus: ",lczetaminus[[i]],")"],{i,len}];];
Do[lcaplus[[i]]=lclambdaplus[[i]] Cos[lczetaplus[[i]]];
lcbplus[[i]]=-lclambdaplus[[i]] Sin[lczetaplus[[i]]];
lccplus[[i]]=lclambdaplus[[i]] Sin[lczetaplus[[i]]];
lcdplus[[i]]=lclambdaplus[[i]] Cos[lczetaplus[[i]]];
lcaminus[[i]]=lclambdaminus[[i]] Cos[lczetaminus[[i]]];
lcbminus[[i]]=-lclambdaminus[[i]] Sin[lczetaminus[[i]]];
lccminus[[i]]=-lclambdaminus[[i]] Sin[lczetaminus[[i]]];
lcdminus[[i]]=-lclambdaminus[[i]] Cos[lczetaminus[[i]]],{i,len}];

lcL=First@Last@Reap[
Do[
Sow[Sqrt[lclambdaplus[[j]]^2+lclambdaminus[[j]]^2+
2*lclambdaplus[[j]] lclambdaminus[[j]] Cos[lczetaplus[[j]]-lczetaminus[[j]]-2 lczetaplus[[1]]]]],{j,len}]
];

If[True,Print["\nLn Scalar:\n ========================== "];
Do[Print["LC-EFA mode: ",i," \tLn = ",lcL[[i]]],{i,len}]];

{lcL,Transpose[{lclambdaplus,lclambdaminus,lczetaplus,lczetaminus}]} 
];


Options[EllipticalFourierDescriptor]={"normalize"->False,"sizeInvariance"->True,"order"->10,"ptsToplot"->300,
"locus"->{0,0},"lobeContribution" -> False};
EllipticalFourierDescriptor[cont:_Image|{{__Integer}..},OptionsPattern[]]:=Module[{contour,coeffL,struct},
contour=generateContour[cont]; (* generate contour *)
coeffL=ellipticalCoefficients[contour,OptionValue["order"]]; (* compute elliptical Coefficients*)
If[OptionValue["normalize"],
coeffL=normalizeEFD[coeffL,OptionValue["sizeInvariance"]]
]; (* to normalize ? *)
struct=plotEFD[contour,coeffL]; (* plot EFD *)
If[OptionValue["lobeContribution"],Join[#,lobeContribution[contour,coeffL]],#]&@{coeffL,struct}
];


End[];


(* ::Section:: *)
(*End Package*)


EndPackage[];

(* ::Package:: *)

harmonicAnalysis[a_,b_,c_,d_]:=Block[{p=a,q=b,r=c,s=d,coeff,\[Theta],rotationMat,starMat,\[Psi],L,W,\[Phi]},
\[Theta] = 0.5ArcTan[(p^2  +r^2 -q^2-s^2),2(p q + r s)];
coeff = {p,q,r,s};
rotationMat = {{Cos[\[Theta]],-Sin[\[Theta]]},{Sin[\[Theta]],Cos[\[Theta]]}};
starMat=Partition[coeff,2,2].rotationMat;
\[Psi]=ArcTan[starMat[[1,1]],starMat[[2,1]]]*(180 /Pi);
\[Phi]=ArcTan[starMat[[1,2]],starMat[[2,2]]]*(180 /Pi);
L = Norm[{starMat[[1,1]],starMat[[2,1]]}];
W= Norm[{starMat[[1,2]],starMat[[2,2]]}];
{L,W,\[Psi],\[Phi]}
 ]


ClearAll@EllipticalFourier;
Options[EllipticalFourier] = {normalize-> True,sizeInvariance-> True,order-> 10,ptsToplot-> 300,locus-> {0.,0.},areathreshold-> 2};
EllipticalFourier[img_,OptionsPattern[]]:=Module[{dXY,dt,t,T,\[Phi],contourpositions, contour,contourtour,contourregion,xi,A0,
\[Delta],C0,dCval,coeffL,normalizeEFD,ellipticalCoefficients,plotEFD,dcCoefficients,dcVal},
contourpositions = PixelValuePositions[MorphologicalPerimeter[img],1];
contourtour=Last@FindShortestTour[contourpositions];
contour = Part[contourpositions,contourtour];
contourregion=BoundaryMeshRegion[contourpositions,Line@contourtour];
dXY = Differences[contour];
dt=(x\[Function]N@Sqrt[x])[Total/@(dXY^2)];
t = Prepend[Accumulate[dt],0];
T= Last@t;
\[Phi] = 2 \[Pi] t/T;

normalizeEFD[coeff_]:= With[
{\[Theta]=0.5 ArcTan[(coeff[[1,1]]^2 -coeff[[1,2]]^2 + coeff[[1,3]]^2 -coeff[[1,4]]^2),
2 (coeff[[1,1]] coeff[[1,2]] +coeff[[1,3]] coeff[[1,4]])]},
Block[{coeffList = coeff,\[Psi],\[Psi]r,a,c,\[Theta]r},
\[Theta]r[i_Integer]:={{Cos[i \[Theta]],-Sin[i \[Theta]]},{Sin[i \[Theta]],Cos[i \[Theta]]}};
coeffList=Partition[#,2,2]&/@coeffList;
coeffList=MapIndexed[(#1.\[Theta]r[First@#2])&,coeffList];
a=First@Flatten@First@coeffList;
c=(Flatten@First@coeffList)[[3]];
\[Psi] = ArcTan[a,c];
\[Psi]r= {{Cos[\[Psi]],Sin[\[Psi]]},{-Sin[\[Psi]],Cos[\[Psi]]}};
coeffList=MapIndexed[Flatten[\[Psi]r.#]&,coeffList];
If[OptionValue[sizeInvariance],a= First@First@coeffList;coeffList/=a,coeffList]]
];

ellipticalCoefficients[order_]:= Module[{coeffs,func},
func[list:{{__}..},n_Integer]:=Block[{listtemp =list, const = T/(2 (n^2)( \[Pi]^2)),phiN= \[Phi] *n,
dCosPhiN, dSinPhiN,an,bn,cn,dn},
dCosPhiN = Cos[Rest@phiN]-Cos[Most@phiN];
dSinPhiN = Sin[Rest@phiN]-Sin[Most@phiN];
an = const*Total@(dXY[[All,1]]/dt *dCosPhiN);
bn = const*Total@(dXY[[All,1]]/dt*dSinPhiN);
cn = const*Total@(dXY[[All,2]]/dt*dCosPhiN);
dn = const*Total@(dXY[[All,2]]/dt*dSinPhiN);
listtemp[[n]]= {an,bn,cn,dn};
listtemp
];
coeffs = ConstantArray[0,{order,4}];
coeffL= If[OptionValue@normalize,normalizeEFD[#],#]&@Fold[func,coeffs,Range@order]
];

dcCoefficients:=( 
xi =  Accumulate[Part[dXY,All,1]]- (dXY[[All,1]]/dt)*Rest@t;
A0= (1/T)*Total@(dXY[[All,1]]/(2 dt) * Differences[t^2] +  xi*dt);
\[Delta] = Accumulate[Part[dXY,All,2]]- (dXY[[All,2]]/dt)*Rest@t;
C0= (1/T)*Total@(dXY[[All,2]]/(2 dt) * Differences[t^2] +  \[Delta]*dt);
{First@First@contour + A0,Last@First@contour+C0}
);

plotEFD[coeff_]:= Last@Reap@Block[{n = Length@coeff,nhalf,nrow=2,ti,xt,yt,func,\[Theta],
i=OptionValue[ptsToplot],j,regionpts,shortesttour,regiondiff,regiondiffmetric},
nhalf = Ceiling[n/2];
ti = Array[#&,i,{0.,1}];
If[Head@OptionValue[locus]===String,dcVal=dcCoefficients,dcVal=OptionValue@locus];
xt = ConstantArray[dcVal[[1]],i];
yt = ConstantArray[dcVal[[2]],i];
regiondiffmetric =100;
j=1;
While[regiondiffmetric>OptionValue[areathreshold],
xt+= coeff[[j,1]] Cos[2 (j) \[Pi] ti] + coeff[[j,2]] Sin[2 (j) \[Pi] ti];
yt += coeff[[j,3]]Cos[2(j) \[Pi] ti] + coeff[[j,4]] Sin[2 (j)\[Pi] ti];

regionpts=Most@Thread@{xt,yt};
shortesttour = Line@Last@FindShortestTour@regionpts;
regiondiff=BoundaryMeshRegion[regionpts,shortesttour]~RegionDifference~contourregion;
regiondiffmetric = Area[regiondiff]/Area[contourregion] * 100;
Sow@{harmonicAnalysis[Sequence@@coeff[[j]]],Graphics[{{Thick,XYZColor[0,0,1,0.6],Line@Thread[{xt,yt}]},
{Thick,Dashed,XYZColor[1,0,0,0.6],Line@contour},Inset[Style[#,{Black,FontSize->12}]&@Text["harmonic: " <> ToString@j]]}]};
j++
] ;
];
ellipticalCoefficients[OptionValue[order]];
plotEFD[coeffL]
]

(*Begin Package*)
BeginPackage["FourierEllipticDescriptors`"];

(*Main Function*)
EllipticFourierDescriptor::usage="
EllipticFourierDescriptor[img_,opts_]
yields the elliptical descriptors for a given binarized image.
Arguments:
img_ = binarized image.
opts_ = different options for which the default values are set as follows:
normalize -> True, sizeInvariance -> True, order -> 10, ptsToplot ->300, locus -> {0.,0.}
";

Options[EllipticFourierDescriptor] = {"normalize"-> True,"sizeInvariance"-> True,"order"-> 10,"ptsToplot"-> 300,"locus"-> {0.,0.}};
EllipticFourierDescriptor[img_,OptionsPattern[]]:=Module[{dXY,dt,t,T,\[Phi],contourpositions, contour,xi,A0,\[Delta],C0,
coeffL,dcCoefficients,dcVal,normalizeEFD,ellipticalCoefficients,plotEFD},
contourpositions = PixelValuePositions[MorphologicalPerimeter[img],1];
contour = Part[contourpositions,Last@FindShortestTour[contourpositions]];
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
If[OptionValue["sizeInvariance"],a= First@First@coeffList;coeffList/=a,coeffList]]
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
coeffL= If[OptionValue@"normalize",normalizeEFD[#],#]&@Fold[func,coeffs,Range@order]
];

dcCoefficients:= ( 
xi =  Accumulate[Part[dXY,All,1]]- (dXY[[All,1]]/dt)*Rest@t;
A0= (1/T)*Total@(dXY[[All,1]]/(2 dt) * Differences[t^2] +  xi*dt);
\[Delta] = Accumulate[Part[dXY,All,2]]- (dXY[[All,2]]/dt)*Rest@t;
C0= (1/T)*Total@(dXY[[All,2]]/(2 dt) * Differences[t^2] +  \[Delta]*dt);
{First@First@contour + A0,Last@First@contour+C0}
);

plotEFD[coeff_]:= Last@Reap@Block[{ti,xt,yt,i=OptionValue["ptsToplot"]},
ti = Array[#&,i,{0.,1}];
If[Head@OptionValue["locus"]===String, dcVal=dcCoefficients, dcVal=OptionValue@"locus"];
xt = ConstantArray[dcVal[[1]],i];
yt = ConstantArray[dcVal[[2]],i];
Do[
xt += coeff[[j,1]] Cos[2 (j) \[Pi] ti] + coeff[[j,2]] Sin[2 (j) \[Pi] ti];
yt += coeff[[j,3]]Cos[2 (j) \[Pi] ti] + coeff[[j,4]] Sin[2 (j)\[Pi] ti];
Sow@Graphics[{{Thick,XYZColor[0,0,1,0.8],Line@contour},{Thick,XYZColor[1,0,0,0.6],Line@Thread[{xt,yt}]},
Inset[Style[#,Bold,Darker@Green,FontSize->12]&@Text["Harmonic: " <> ToString@j]]}];
,{j,1,Length@coeff}];
];

ellipticalCoefficients[OptionValue["order"]];
plotEFD[coeffL]
];

(*Private Content*)
Begin["`Private`"];
End[];

(*End Package*)
EndPackage[];

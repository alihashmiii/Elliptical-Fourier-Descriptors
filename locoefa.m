(* ::Package:: *)

(* ::Section:: *)
(*Begin Package*)


BeginPackage["LobeContributionEFA`"]


(* ::Section:: *)
(*Main Function*)


Options[EllipticalFourier] = {normalize -> False, 
   sizeInvariance -> True, order -> 10, ptsToplot -> 300, 
   locus -> "dcCoefficients"};
EllipticalFourier[cont : _Image | {{__Integer} ..}, 
  OptionsPattern[]] := Module[{dXY, dt, t, T, \[Phi], contourpositions, contour, contourtour, contourregio, xi, A0, \[Delta], C0, dCval, coeffL,
   normalizeEFD, ellipticalCoefficients, plotEFD, struct},
  If[Head@cont === Image, contourpositions = DeleteDuplicates@PixelValuePositions[MorphologicalPerimeter[cont], 1];,
   contourpositions = DeleteDuplicates@cont;
   ];
  contourtour = Last@FindShortestTour[contourpositions];
  contour = Part[contourpositions, contourtour];
  dXY = Differences[contour];
  dt = (x \[Function] N@Sqrt[x])[Total /@ (dXY^2)];
  t = Prepend[Accumulate[dt], 0];
  T = Last@t; \[Phi] = 2 \[Pi] t/T;
  
  normalizeEFD[coeff_] := With[{\[Theta] = 0.5 ArcTan[(coeff[[1, 1]]^2 - coeff[[1, 2]]^2 + 
          coeff[[1, 3]]^2 - coeff[[1, 4]]^2), 2 (coeff[[1, 1]] coeff[[1, 2]] + coeff[[1, 3]] coeff[[1, 4]])]},
    Block[{coeffList = coeff, \[Psi], \[Psi]r, a, c, \[Theta]r},
     \[Theta]r[i_Integer] := {{Cos[i \[Theta]], -Sin[i \[Theta]]}, {Sin[i \[Theta]], Cos[i \[Theta]]}};
     coeffList = Partition[#, 2, 2] & /@ coeffList;
     coeffList = MapIndexed[(#1.\[Theta]r[First@#2]) &, coeffList];
     a = First@Flatten@First@coeffList;
     c = (Flatten@First@coeffList)[[3]];
     \[Psi] = ArcTan[a, c];
     \[Psi]r = {{Cos[\[Psi]], Sin[\[Psi]]}, {-Sin[\[Psi]], Cos[\[Psi]]}};
     coeffList = MapIndexed[Flatten[\[Psi]r.#] &, coeffList];
     If[OptionValue[sizeInvariance], a = First@First@coeffList; 
      coeffList /= a, coeffList]]
    ];
  
  ellipticalCoefficients[order_] := Module[{coeffs, func},
    func[list : {{__} ..}, n_Integer] := Block[{listtemp = list, const = T/(2 (n^2) ( \[Pi]^2)), 
       phiN = \[Phi] *n, dCosPhiN, dSinPhiN, an, bn, cn, dn},
      dCosPhiN = Cos[Rest@phiN] - Cos[Most@phiN];
      dSinPhiN = Sin[Rest@phiN] - Sin[Most@phiN];
      an = const*Total@(dXY[[All, 1]]/dt *dCosPhiN);
      bn = const*Total@(dXY[[All, 1]]/dt*dSinPhiN);
      cn = const*Total@(dXY[[All, 2]]/dt*dCosPhiN);
      dn = const*Total@(dXY[[All, 2]]/dt*dSinPhiN);
      listtemp[[n]] = {an, bn, cn, dn};
      listtemp
      ];
    coeffs = ConstantArray[0, {order, 4}];
    coeffL = If[OptionValue@normalize, normalizeEFD[#], #]&@Fold[func, coeffs, Range@order]
    ];
  
  dcCoefficients := (
    xi =  Accumulate[Part[dXY, All, 1]] - (dXY[[All, 1]]/dt)*Rest@t;
    A0 = (1/T)* Total@(dXY[[All, 1]]/(2 dt) * Differences[t^2] +  xi*dt);
    \[Delta] = Accumulate[Part[dXY, All, 2]] - (dXY[[All, 2]]/dt)*Rest@t;
    C0 = (1/T)* Total@(dXY[[All, 2]]/(2 dt) * Differences[t^2] +  \[Delta]*dt);
    {First@First@contour + A0, Last@First@contour + C0}
    );
  
  plotEFD[coeff_] := First@Last@Reap@Block[{n = Length@coeff, nhalf, nrow = 2, ti, xt, yt, 
        func, \[Theta], i = OptionValue[ptsToplot]},
       nhalf = Ceiling[n/2];
       ti = Array[# &, i, {0., 1}];
       If[Head@OptionValue[locus] === String, dcVal = dcCoefficients, dcVal = OptionValue@locus];
       xt = ConstantArray[dcVal[[1]], i];
       yt = ConstantArray[dcVal[[2]], i];
       Do[
        xt += coeff[[j, 1]] Cos[2 (j) \[Pi] ti] + coeff[[j, 2]] Sin[2 (j) \[Pi] ti];
        yt += coeff[[j, 3]] Cos[2 (j) \[Pi] ti] + coeff[[j, 4]] Sin[2 (j) \[Pi] ti];
        Sow[{coeff[[j]], Graphics[{{Thick, XYZColor[0, 0, 1, 0.6], Line@Thread[{xt, yt}]}, {Thick, Dashed, Black, 
             Line@contour}, Inset[Style[#, {XYZColor[1, 0, 0, 0.6], Bold, FontSize -> 12}]&@Text["harmonic: "
         <> ToString@j]]}]}]
         ,{j,OptionValue[order]}]
       ];
  ellipticalCoefficients[OptionValue[order]];
  struct = plotEFD[coeffL][[All, -1]];
  {struct, LobeContributionEFA[contour, coeffL, False]}
  ]


LobeContributionEFA[contour : {{__Integer} ..}, modes_, debug_: False] := 
 Module[{starmat = {}, lcL = {}, lambdaminus, lambdaminusvec = {}, lambdaplus, lambdaplusvec = {}, zetaplus, zetaplusvec = {}, 
   zetaminus, zetaminusvec = {}, lczetaplus, lclambdaplus, lczetaminus, lclambdaminus, lcaplus, lcbplus, lccplus, lcdplus, 
   lcaminus, lcbminus, lccminus, lcdminus, maxerror = 1 10^-13, tau, alphaprime, gammaprime, \[Rho], alphastar, betastar,
   gammastar, deltastar, r, a, b, c, d, \[Phi], aprime, bprime, cprime, dprime, \[Theta], lambda1, lambda12, lambda21, lambda2,
   lclambdaplusZeroM, lczetaplusZeroM, len = Length@modes},
  
  lczetaplus = lclambdaplus = lczetaminus = lclambdaminus = ConstantArray[0, len];
  lcaplus = lcbplus = lccplus = lcdplus = lcaminus = lcbminus = lccminus = lcdminus =  ConstantArray[0, len];
  
  tau = 0.5 ArcTan[modes[[1, 1]]^2 + modes[[1, 3]]^2 - modes[[1, 2]]^2 - modes[[1, -1]]^2, 
     2.0*(modes[[1, 1]] modes[[1, 2]] + modes[[1, 3]] modes[[1, -1]])];
  alphaprime = modes[[1, 1]] Cos[tau] + modes[[1, 2]] Sin[tau];
  gammaprime = modes[[1, 3]] Cos[tau] + modes[[1, -1]] Sin[tau];
  \[Rho] = ArcTan[alphaprime, gammaprime];
  If[\[Rho] < 0, tau += N[\[Pi]]];
  Do[
   alphastar = modes[[k, 1]] Cos[k tau] + modes[[k, 2]] Sin[k tau];
   betastar = - modes[[k, 1]] Sin[k tau] + modes[[k, 2]] Cos[k tau];
   gammastar = modes[[k, 3]] Cos[k tau] + modes[[k, -1]] Sin[k tau];
   deltastar = - modes[[k, 3]] Sin[k tau] + modes[[k, -1]] Cos[k tau];
   AppendTo[starmat, {alphastar, betastar, gammastar, deltastar}];
   , {k, len}];
  r = starmat[[1, 1]] starmat[[1, -1]] - starmat[[1, 2]] starmat[[1, -2]];
  If[r < 0, starmat[[All, {2, -1}]] *= -1];
  {a, b, c, d} = starmat\[Transpose];
  
  If[debug,
   Print["modified EFA coefficients\n =============== "];
   Do[Print["mode:", i, {a[[i]], b[[i]], c[[i]], d[[i]]}], {i, len}];
   ];
  If[debug, Print["lambda matrices:\n =============== "]];
  Do[
   \[Phi] = 0.5 ArcTan[a[[i]]^2 + c[[i]]^2 - b[[i]]^2 - d[[i]]^2, 2.0 (a[[i]] b[[i]] + c[[i]] d[[i]])];
   aprime = a[[i]] Cos[\[Phi]] + b[[i]] Sin[\[Phi]];
   bprime = - a[[i]] Sin[\[Phi]] + b[[i]] Cos[\[Phi]];
   cprime = c[[i]] Cos[\[Phi]] + d[[i]] Sin[\[Phi]];
   dprime = -c[[i]] Sin[\[Phi]] + d[[i]] Cos[\[Phi]];
   \[Theta] = ArcTan[aprime, cprime];
   lambda1 = aprime Cos[\[Theta]] + cprime Sin[\[Theta]];
   lambda12 = bprime Cos[\[Theta]] + dprime Sin[\[Theta]];
   lambda21 = - aprime Sin[\[Theta]] + cprime Cos[\[Theta]];
   lambda2 = - bprime Sin[\[Theta]] + dprime Cos[\[Theta]];
   If[Or[debug, Abs[lambda12] > maxerror, Abs[lambda21] > maxerror, lambda1 < 0, lambda1 < Abs@lambda2],
    Which[Abs[lambda12] > maxerror || Abs[lambda21] > maxerror,
      Print["off-diagonal lambda matrix unequal to zero:\n"],
      lambda1 < 0, Print["Warning: lambda1 negative\n"],
      lambda1 < Abs[lambda2], 
      Print["Warning: lambda 1 < |lambda2|:\n"];
      Print["mode:\n", "(lambda1,lambda12): ", {lambda1, lambda12}, "\n", "(lambda21,lambda2): ", {lambda1, lambda12}, "\n"];
        ];
    ];
   lambdaplus = (lambda1 + lambda2)/2.0;
   lambdaminus = (lambda1 - lambda2)/2.0;
   zetaplus = \[Theta] - \[Phi];
   zetaminus = -\[Theta] - \[Phi];
   AppendTo[zetaplusvec, zetaplus];
   AppendTo[lambdaplusvec, lambdaplus];
   AppendTo[lambdaminusvec, lambdaminus];
   AppendTo[zetaminusvec, zetaminus];, {i, Length@starmat}];
  
  lclambdaplusZeroM = lambdaplusvec[[2]];
  lczetaplusZeroM = zetaplusvec[[2]];
  
  lclambdaplus[[1]] = lambdaplusvec[[1]];
  lczetaplus[[1]] = zetaplusvec[[1]];
  Do[
   lclambdaplus[[j]] = lambdaplusvec[[j + 1]];
   lczetaplus[[j]] = zetaplusvec[[j + 1]]
   , {j, 2, len - 1}];
  Do[
   lclambdaminus[[j]] = lambdaminusvec[[j - 1]];
   lczetaminus[[j]] = zetaminusvec[[j - 1]]
   , {j, 2, len}];
  If[True, Print["Ln quadruplets:\n ======================"];
   Print["lc-EFA mode 0: ", {lclambdaplusZeroM, lczetaplusZeroM}];
   Do[Print["lc-EFA mode ", i, "\n", "(lambdaplus: ", lclambdaplus[[i]], " lambdaminus: ", lclambdaminus[[i]], ")\n", 
     "(zetaplus: ", lczetaplus[[i]], " zetaminus: ", lczetaminus[[i]], ")"]
    , {i, len}];
   ];
  Do[
   lcaplus[[i]] = lclambdaplus[[i]] Cos[lczetaplus[[i]]];
   lcbplus[[i]] = -lclambdaplus[[i]] Sin[lczetaplus[[i]]]; 
   lccplus[[i]] = lclambdaplus[[i]] Sin[lczetaplus[[i]]]; 
   lcdplus[[i]] = lclambdaplus[[i]] Cos[lczetaplus[[i]]]; 
   lcaminus[[i]] = lclambdaminus[[i]] Cos[lczetaminus[[i]]]; 
   lcbminus[[i]] = -lclambdaminus[[i]] Sin[lczetaminus[[i]]]; 
   lccminus[[i]] = -lclambdaminus[[i]] Sin[lczetaminus[[i]]]; 
   lcdminus[[i]] = -lclambdaminus[[i]] Cos[lczetaminus[[i]]];
   , {i, len}];
  If[debug, Print["LC coefficients:\n======================"];
   Do[Print["\nmode: ", i, "A+\n", {lcaplus[[i]], lcbplus[[i]]}, 
     "\n", {lccplus[[i]], lcdplus[[i]]}];
    Print["\nmode: ", i, "A-\n", {lcaminus[[i]], lcbminus[[i]]}, "\n", {lccminus[[i]], lcdminus[[i]]}];
    , {i, len}];
   ];
  Do[
  AppendTo[lcL,
   Sqrt[lclambdaplus[[j]]^2 + lclambdaminus[[j]]^2 + 2*lclambdaplus[[j]] lclambdaminus[[j]] Cos[
        lczetaplus[[j]] - lczetaminus[[j]] - 2 lczetaplus[[1]]]]]
        , {j, len}];
  
  If[True, Print["\nLn Scalar:\n ========================== "];
   Do[Print["LC-EFA mode: ", i, " \tLn = ", lcL[[i]]], {i, len}]];
  {lcL, Transpose[{lclambdaplus, lclambdaminus, lczetaplus, 
     lczetaminus}]}
  ]


(* ::Section:: *)
(*Private Content*)


Begin["`Private`"];
End[];


(* ::Section:: *)
(*End Package*)


EndPackage[];

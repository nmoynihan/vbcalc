(* ::Package:: *)

BeginPackage["VBCalc`"];

version = 0.1;

VB::usage = "VB[metric, coordinates, assumptions] returns a vierbein of a given metric. The metric input must be square and symmetric.";
BVector::usage = "BVector[vierbein, assumptions] calculates the vector \!\(\*SuperscriptBox[\(B\), \(d\)]\)=\!\(\*SuperscriptBox[\(\[Epsilon]\), \(abcd\)]\)\!\(\*SubscriptBox[\(\[Omega]\), \(bca\)]\), where \!\(\*SubscriptBox[\(\[Omega]\), \(bca\)]\) is the spin connection.";
CompareVB::usage = "CompareVB[metric, vierbein1, vierbein2, assumptions] calculates the Lorentz transformation that connects two related vierbeins, of the form \!\(\*SubscriptBox[SuperscriptBox[\(e\), \(a'\)], \(\[Mu]\)]\)=\!\(\*SubscriptBox[SuperscriptBox[\(\[CapitalLambda]\), \(a'\)], \(a\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(e\), \(a\)], \(\[Mu]\)]\). If the two vierbeins are not related by a Lorentz transform, (i.e. det(\[CapitalLambda]) is not +1 or -1, or \!\(\*SuperscriptBox[\(\[CapitalLambda]\), \(T\)]\)\[Eta]\[CapitalLambda] != \[Eta]) then an error is printed.";
CalcET::usage = "CalcET[metric, coordinates, print = 1 or 0, assumptions] returns a vector containing the Affine Connection, Riemann tensor, Ricci Tensor, Scalar Curvature and Einstein Tensor. Setting print = 1 also prints the non-zero components of each of these.";

Begin["`Private`"];

VB[metric_, cd_, as_] :=
Module[{i,k,m,B,n,L},
If[Dimensions[metric][[1]] === Dimensions[metric][[2]] ,Null,Print["The metric input is not square. Only square matrices can be processed"];Abort[]]; (* Check if square *)
If[SymmetricMatrixQ[metric],Null,Print["The metric input is not symmetric, please input a square, symmetric matrix."]; Abort[]]; (* Check if Symmetric *)
B = metric;
B[[1,1]] = -B[[1,1]];
n = Dimensions[B][[1]];
L = Table[0,{n},{n}];
For[k=1,k<= n, k++,
L[[k,k]]= Sqrt[B [[k,k]]-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(m = 1\), \(k - 1\)]\(
\*SuperscriptBox[\((L[\([k, m]\)])\), \(2\)]*If[m == 1, \(-1\), 1]\)\)];
(* If[k==1,Print["Debug: ", (B [[i,k]]-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(m = 1\), \(k - 1\)]\(L[\([i, m]\)] L[\([k, m]\)]\)\))/L[[k,k]] ]; Print["Debug k: ", k ],1]; *)
For[i=k+1,i<= n, i++,
(* If[k==1,Print["Debug: ", (B [[i,k]]-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(m = 1\), \(k - 1\)]\(L[\([i, m]\)] L[\([k, m]\)]*If[m == 1, \(-1\), 1]\)\))/L[[k,k]] ]; Print["Debug k: ", k ],1]; *)
L[[i,k]]=If[k==1,-1,1] (B [[i,k]]-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(m = 1\), \(k - 1\)]\(L[\([i, m]\)] L[\([k, m]\)]*If[m == 1, \(-1\), 1]\)\))/L[[k,k]]
]
];
vb = Transpose[L]//Simplify;
inversevb = Inverse[vb];
eta = ({
 {-1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}
});
coord = cd;
result = L.eta.vb- metric//FullSimplify;
newmetric = L.eta.vb//FullSimplify;
(* If[result ===Table[0,{n},{n}],Print["This vierbein passes the test and looks to be correct:"],Print["This vierbein did not pass the test and may be incorrect. (DEBUG:) Result = ", MatrixForm[result]]]; *)
dcoord="d"<>ToString[#]&/@coord//ToExpression/@#& ;
le = dcoord.newmetric.dcoord//FullSimplify;
Print["From this vierbein, we recover the line Element = ", FactorTerms[le,dcoord]];
If[latex ===1,Print["Line Element: ", ToString[FactorTerms[le,dcoord], TeXForm]],Null];
FullSimplify[vb, Assumptions-> as]
];

End[];


Begin["`Private`"];

<<<<<<< HEAD
BVector[vb_, coords_, as_] := Module[{i,inversevb},
inversevb = Inverse[vb];
B =Table[Sum[LeviCivitaTensor[4][[a,b,c,i]]eta[[b,k]]vb[[k,\[Lambda]]]D[inversevb[[\[Lambda],c]],coords[[a]]],{a,4},{b,4},{c,4},{\[Lambda],4},{k,4}],{i,1,4}]; (* Calculate the B0 curvature term *)
=======
BVector[vb_, as_] := Module[{i},
B =Table[Sum[LeviCivitaTensor[4][[a,b,c,i]]eta[[b,k]]vb[[k,\[Lambda]]]D[inversevb[[\[Lambda],c]],coord[[a]]],{a,4},{b,4},{c,4},{\[Lambda],4},{k,4}],{i,1,4}]; (* Calculate the B0 curvature term *)
>>>>>>> origin/master
FullSimplify[B, Assumptions-> as]
]

End[];


Begin["`Private`"];

CompareVB[metric_, vb1_, vb2_, as_] :=
Module[{i,k,m},
L = vb1.metric.Transpose[vb2];
L1 = vb1.Inverse[vb2];
res = L + L1;
eta = ({
 {-1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}
});
L2 = FullSimplify[L1, Assumptions-> as];
(* Print["Lorentz Transformation = ", L2//MatrixForm//FullSimplify];
Print["Result: ", res//FullSimplify//MatrixForm];
Print["This doesn't seem to be working...."]; *)
(* This Determinant should be +/- 1 *)
dt = FullSimplify[Det[L2], Assumptions-> as];
Print["Determinant of LT: ", dt];
(* This actually produces an inverse LT, which is what L2 is *)
check = FullSimplify[Transpose[Inverse[L2]].eta.Inverse[L2], Assumptions->as];
Print["\!\(\*SuperscriptBox[\(\[CapitalLambda]\), \(T\)]\)\[Eta]\[CapitalLambda] = ", check//MatrixForm];
If [check=== eta \[And]  dt ===1 \[Or] dt === -1, Print["This looks like a valid Lorentz Transformation"], Print["Can't compute a valid Lorentz Transformation, vierbeins may not be equivalent."]; Abort[]];
L2
]

End[];



Begin["`Private`"];

<<<<<<< HEAD
CalcET[metric_, coord_, printet_, as_] := Module[{affine, riemann, inversemetric, et, scalar,listaffine,listriemann,n,i,j,k,R,ricci,listricci,einstein,listeinstein,output},
inversemetric = Inverse[metric];
n=4;
et = {0,0,0,0,0};
=======
CalcET[metric_, coord_, printet_, as_] := Module[{affine, riemann, inversemetric, et, scalar},
inversemetric = Inverse[metric];
>>>>>>> origin/master
affine = Simplify[Table[(1/2)*Sum[(inversemetric[[i,s]])*(D[metric[[s,j]],coord[[k]] ]+D[metric[[s,k]],coord[[j]] ]-D[metric[[j,k]],coord[[s]] ]),{s,1,n}],{i,1,n},{j,1,n},{k,1,n}] ];
et[[1]] = affine;
listaffine := Table[If[UnsameQ[affine[[i,j,k]],0],{ToString[\[CapitalGamma][i,j,k]],affine[[i,j,k]]}] ,{i,1,n},{j,1,n},{k,1,j}];
If[printet === 1,Print["Affine Connection:"]; Print[TableForm[Partition[DeleteCases[Flatten[listaffine],Null],2],TableSpacing->{2,2}]];,Null];
riemann := Simplify[Table[D[affine[[i,j,l]],coord[[k]] ]-D[affine[[i,j,k]],coord[[l]] ]+Sum[affine[[s,j,l]] affine[[i,k,s]]-affine[[s,j,k]] affine[[i,l,s]],{s,1,n}],{i,1,n},{j,1,n},{k,1,n},{l,1,n}] ];
et[[2]] = riemann;
listriemann:=Table[If[UnsameQ[riemann[[i,j,k,l]],0],{ToString[R[i,j,k,l]],riemann[[i,j,k,l]]}] ,{i,1,n},{j,1,n},{k,1,n},{l,1,k-1}];
If[printet === 1, Print["Riemann Tensor:"]; Print[TableForm[Partition[DeleteCases[Flatten[listriemann],Null],2],TableSpacing->{2,2}]];,Null];
ricci:=ricci=Simplify[Table[Sum[riemann[[i,j,i,l]],{i,1,n}],{j,1,n},{l,1,n}] ];
et[[3]] = ricci;
listricci:=Table[If[UnsameQ[ricci[[j,l]],0],{ToString[R[j,l]],ricci[[j,l]]}] ,{j,1,n},{l,1,j}];
If[printet === 1, Print["Ricci Tensor:"];Print[TableForm[Partition[DeleteCases[Flatten[listricci],Null],2],TableSpacing->{2,2}]];,Null];
scalar=Simplify[Sum[inversemetric[[i,j]]ricci[[i,j]],{i,1,n},{j,1,n}]];
et[[4]] = scalar;
<<<<<<< HEAD
If[printet === 1, Print["Scalar Curvature, R =", scalar];,Null];
=======
Ifprintet === 1, Print["Scalar Curvature, R =", scalar];,Null];
>>>>>>> origin/master
einstein:=einstein=Simplify[ricci-(1/2)scalar*metric];
et[[5]] = einstein;
listeinstein:=Table[If[UnsameQ[einstein[[j,l]],0],{ToString[G[j,l]],einstein[[j,l]]}] ,{j,1,n},{l,1,j}];
If[printet === 1, Print["Einstein Tensor:"]; Print[TableForm[Partition[DeleteCases[Flatten[listeinstein],Null],2],TableSpacing->{2,2}]];,Null];
output = et
]

End[];

Begin["`Private`"];

Print["VBCalc ", version, " Loaded and ready."];

End[];

EndPackage[];

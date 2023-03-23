(* ::Package:: *)

BeginPackage["GaussianPiecewiseIntegrals`"];


GaussianIntegral::usage = "Integrate from -\[Infinity] to \[Infinity], speeding up Gaussian integrals. Defaults to Integrate when it cannot identify a Gaussian integral. Also allows for integration of Hermite polynomials multiplied by Erf[alpha x+beta].";


Begin["`Private`"]


IsElExp[el_] :=
	(Head[el] === Power && Length[el] > 1 && el[[1]] === E);

IsExp[list_] := IsElExp /@ list;
	
MakeFactorList[expr_]:=Module[{factorlist},
	If[expr===0||expr===0.0,Return[0]];
	factorlist = Replace[expr, Times -> List, 1, Heads -> True];
	If[Not[Head[factorlist] === List],
		factorlist = {factorlist}
	];
	factorlist
];
	
GetCoeffs[factorlist_,intvar_]:=Module[{mask,exp,prefactor,coeffs,a,b,c,n,const},
	mask = IsExp[factorlist];
		If[Not[AnyTrue[mask, (# == True)&]],
			(
				Message[SingleGaussianIntegral::mathematica, expr];
				Return[Integrate[expr, {intvar, -Infinity, Infinity}], Module]
			)
		];
		exp = Expand[(Times @@ Pick[factorlist, mask])[[2]]];
		prefactor = Times @@ Pick[factorlist, Not /@ mask];
		coeffs = CoefficientList[exp, intvar];
		If[Length[coeffs] > 3,
			(
				Message[SingleGaussianIntegral::mathematica, expr];
				Return[Integrate[expr, {intvar, -Infinity, Infinity}], Module]
			)
		];
		{c, b, a} = coeffs;
		n = Exponent[prefactor, intvar];
		const = Coefficient[prefactor, intvar, n];
		{const,{a,b,c,n}}
		];
GaussianIntegralResult[a_, b_, c_, n_] :=
	1 / 2 (-a) ^ (-n / 2) E^c (1 / a (-1 + (-1) ^ n) b Gamma[1 + n / 2] 
		Hypergeometric1F1[1 + n / 2, 3 / 2, -(b^2 / (4 a))] + ((1 + (-1) ^ n)
		 Gamma[(1 + n) / 2] Hypergeometric1F1[(1 + n) / 2, 1 / 2, -(b^2 / (4 
		a))]) / Sqrt[-a]);
SingleGaussianIntegral[expr_, intvar_] :=
	Module[{factorlist, coeffs},
		If[expr===0||expr===0.0,Return[0]];
		factorlist=MakeFactorList[expr];
		coeffs=GetCoeffs[factorlist,intvar];
		coeffs[[1]] * GaussianIntegralResult@@coeffs[[2]]
	]

SingleGaussianIntegral::mathematica = "Evaluating `1` using Mathematica. Results may be slow.";
PlainGaussianIntegral[expr_, intvar_, simplify_:False] :=
	Module[{termslist, sum},
		termslist = Replace[Expand[expr], Plus -> List, 1, Heads -> True];
		If[Not[Head[termslist] === List],
			termslist = {termslist}
		];
		sum=0;
		Do[sum+=SingleGaussianIntegral[termslist[[i]],intvar];,{i, 1, Length
	[termslist]}];
		
		(*sum = Total[Parallelize[Map[SingleGaussianIntegral[#, intvar]&, termslist
			
			]]];*)
		If[simplify,
			Simplify[sum]
			,
			sum
		]
	]


IsPiecewise[el_] :=Head[el] === Piecewise;
NegativePiecewiseIntegral[a_,b_,c_,n_]:=1/(2 a) (-1)^n (-a)^(-n/2) E^c (b Gamma[1+n/2] Hypergeometric1F1[1+n/2,3/2,-(b^2/(4 a))]-Sqrt[-a] Gamma[(1+n)/2] Hypergeometric1F1[(1+n)/2,1/2,-(b^2/(4 a))]);
PositivePiecewiseIntegral[a_,b_,c_,n_]:=1/2 (-a)^(-1-n/2) E^c (b Gamma[1+n/2] Hypergeometric1F1[1+n/2,3/2,-(b^2/(4 a))]+Sqrt[-a] Gamma[(1+n)/2] Hypergeometric1F1[(1+n)/2,1/2,-(b^2/(4 a))]);

PiecewiseIntegral[expr_,intvar_]:=Module[{second,negative,positive,valpositive,valnegative},
	positive=expr[[1]];
	second=expr[[2]];
	If[Length[positive]>1,second=positive[[2]]];
	positive=positive[[1]];
	If[positive[[2]]===(intvar>=0)||positive[[2]]===(intvar>0),
		negative=second;
		While[Head[negative]==List,negative=First[negative]];
		positive=positive[[1]];
	,
		negative=positive[[1]];
		positive=second;
		While[Head[positive]==List,positive=First[positive]];
	];
	
	valpositive=If[positive===0||positive===0.0||positive===Indeterminate||positive===Undefined,0,SinglePiecewiseIntegral[positive,intvar,PositivePiecewiseIntegral]];
	valnegative=If[negative===0||negative===0.0||negative===Indeterminate||negative===Undefined,0,SinglePiecewiseIntegral[negative,intvar,NegativePiecewiseIntegral]];
	
	
	
	valpositive+valnegative
];



SinglePiecewiseIntegral[expr_, intvar_, func_] :=
	Module[{termslist,factorlist,coeffs},	
		termslist = Replace[Expand[expr], Plus -> List, 1, Heads -> True];
		If[Not[Head[termslist] === List],
			termslist = {termslist}
		];
		coeffs = GetCoeffs[MakeFactorList[#],intvar]&/@termslist;
		
		Total[#[[1]] * func@@#[[2]]&/@coeffs]
	];

GaussianIntegral[expr_, intvar_, simplify_:False] :=
	Module[{sum, piecewise},
		piecewise=PiecewiseExpand[expr];
		sum = If[IsPiecewise[piecewise],PiecewiseIntegral[piecewise,intvar],PlainGaussianIntegral[piecewise,intvar]];
		If[simplify,
			Simplify[sum]
			,
			sum
		]
	]


End[];
EndPackage[];

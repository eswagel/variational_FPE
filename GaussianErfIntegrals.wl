(* ::Package:: *)

BeginPackage["GaussianErfIntegrals`"];


GaussianIntegral::usage = "Integrate from -\[Infinity] to \[Infinity], speeding up Gaussian integrals. Defaults to Integrate when it cannot identify a Gaussian integral. Also allows for integration of Hermite polynomials multiplied by Erf[alpha x+beta].";


Begin["`Private`"]


W37[n_, alpha_, alpha1_, beta_, beta1_] :=
	Sqrt[\[Pi] / alpha] Exp[beta^2 / (4 alpha)] Erf[(2 alpha beta1 + alpha1 
		beta) / (2 Sqrt[alpha (alpha + alpha1^2)])] Sum[beta ^ (n - 2 l) / (l
		! * (n - 2 l)! alpha ^ (n - l)), {l, 0, Floor[n / 2]}] + 1 / Sqrt[alpha
		 + alpha1^2] Exp[(beta^2 - 4 alpha beta1^2 - 4 alpha1 beta beta1) / (
		4 alpha + 4 alpha1^2)] (Sum[(l! * beta ^ (n + 1 - 2 l)) / ((2 l)! (n 
		+ 1 - 2 l)!) Sum[4 ^ (l + 1 - r) (2 r - 2)! / (r - 1)! Sum[alpha1 ^ (
		2 r - 1 - 2 q) (2 alpha beta1 + alpha1 beta) ^ (2 r - 2 - 2 q) / (q! 
		(2 r - 2 - 2 q)! alpha ^ (n - 1 + r - 2 q) (alpha + alpha1^2) ^ (2 r 
		- 2 - q)), {q, 0, r - 1}], {r, 1, l}], {l, 1, n - Floor[n / 2]}] - Sum[
		beta ^ (n - 2 l) / (l! (n - 2 l)!) Sum[(r - 1)! Sum[alpha1 ^ (2 r - 2
		 q) (2 alpha beta1 + alpha1 beta) ^ (2 r - 1 - 2 q) / (q! (2 r - 1 - 
		2 q)! alpha ^ (n - l + r - 2 q) (alpha + alpha1^2) ^ (2 r - 1 - q)), 
		{q, 0, r - 1}], {r, 1, l}], {l, 1, Floor[n / 2]}])

HermiteErfInt[n_, a_, b_, c_, alpha_, beta_] :=
	n! / 2^n * W37[n, -a, alpha, b, beta] * Exp[c];
	
HermiteErfIntNoB[n_,a_,c_,alpha_]:=(-1+(-1)^n)((-a)^(-n/2)alpha Exp[c] Gamma[1+n/2]Hypergeometric2F1[1/2,(2+n)/2,3/2,alpha^2/a])/(a Sqrt[\[Pi]]);

GaussianIntegralResult[a_, b_, c_, n_] :=
	1 / 2 (-a) ^ (-n / 2) E^c (1 / a (-1 + (-1) ^ n) b Gamma[1 + n / 2] 
		Hypergeometric1F1[1 + n / 2, 3 / 2, -(b^2 / (4 a))] + ((1 + (-1) ^ n)
		 Gamma[(1 + n) / 2] Hypergeometric1F1[(1 + n) / 2, 1 / 2, -(b^2 / (4 
		a))]) / Sqrt[-a]);

IsElExp[el_] :=
	(Head[el] === Power && Length[el] > 1 && el[[1]] === E);

IsExp[list_] :=
	IsElExp /@ list

IsErf[list_, intvar_] :=
	(
			Head[#] === Erf &&
				If[Length[#] > 0,
					Exponent[#[[1]], intvar] > 0
					,
					False
				]
		)& /@ list

SingleIntegral[expr_, intvar_] :=
	Module[{factorlist, exp, coeffs, mask, prefactor, a, b, c, n, const,
		 erfmask, erf, prefactornoerf, erfcoeffs, alpha, beta},
		factorlist = Replace[expr, Times -> List, 1, Heads -> True];
		If[Not[Head[factorlist] === List],
			factorlist = {factorlist}
		];
		mask = IsExp[factorlist];
		If[Not[AnyTrue[mask, (# == True)&]],
			Return[Integrate[expr, {intvar, -Infinity, Infinity}], Module]
		];
		exp = Expand[(Times @@ Pick[factorlist, mask])[[2]]];
		prefactor = Pick[factorlist, Not /@ mask];
		coeffs = CoefficientList[exp, intvar];
		If[Length[coeffs] > 3,
			Message[SingleGaussianIntegral::mathematica, expr];
			Return[Integrate[expr, {intvar, -Infinity, Infinity}]]
		];
		{c, b, a} = coeffs;
		erfmask = IsErf[prefactor, intvar];
		erf = Pick[prefactor, erfmask];
		If[Length[erf] > 1,
			Message[SingleGaussianIntegral::mathematica, expr];
			Return[Integrate[expr, {intvar, -Infinity, Infinity}]];
		];
		prefactornoerf = Times @@ Pick[prefactor, Not /@ erfmask];
		n = Exponent[prefactornoerf, intvar];
		const = Coefficient[prefactornoerf, intvar, n];
		If[Length[erf] > 0,
			(
				erf = erf[[1]][[1]];
				{beta, alpha} = CoefficientList[erf, intvar];
			)
			,
			Return[const * GaussianIntegralResult[a, b, c, n]]
		];
		If[b == 0 && beta == 0,
			Return[const*HermiteErfIntNoB[n,a,c,alpha]];
		];
		const * HermiteErfInt[n, a, b, c, alpha, beta]
	]

SingleIntegral::mathematica = "Evaluating `1` using Mathematica. Results may be slow.";

GaussianIntegral[expr_, intvar_, simplify_:False] :=
	Module[{termslist, sum},
		termslist = Replace[Expand[expr], Plus -> List, 1, Heads -> True];
		If[Not[Head[termslist] === List],
			termslist = {termslist}
		];
		sum = 0;
		Do[
			sum += SingleIntegral[termslist[[i]], intvar];
			,
			{i, Length[termslist]}
		];
		If[simplify,
			Simplify[sum]
			,
			sum
		]
	]


End[];
EndPackage[];

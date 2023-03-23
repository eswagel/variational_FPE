(* ::Package:: *)

(* ::Input::Initialization:: *)
PrependTo[$Path,"C:/Users/ethan/OneDrive/Documents/Wolfram Mathematica"];
<< "GaussianIntegrals`"
IntInf[f_,intvar_:x]:=GaussianIntegral[f,intvar];
calcStats[f_,H_]:=Module[{fl2=f/Sqrt[IntInf[f^2]]},{N[IntInf[fl2*H[fl2]]],N[IntInf[x^2 f]],N[IntInf[x^4 f]]}]
calcStatsNoH[f_]:={N[IntInf[x^2 f]],N[IntInf[x^4 f]]}
reportStats[f_,H_]:=reportStats[calcStats[f,H]]
ToStringMod[x_]:=If[NumericQ[x]&&Im[x]==0&&x!=0&&Re[Log10[x]]<=-5,ToString[ScientificForm[x],TraditionalForm],ToString[x]]
reportStats[stats_]:="\[LeftAngleBracket]x^2\[RightAngleBracket] = "<>ToStringMod[stats[[2]]]<>"\n \[LeftAngleBracket]x^4\[RightAngleBracket] = "<>ToStringMod[stats[[3]]]<>"\n \[LeftAngleBracket]H\[RightAngleBracket] = "<>ToStringMod[stats[[1]]];
reportStatsComparison[stats_,comparison_]:=Module[{del={stats[[1]],stats[[2]]-comparison[[2]],stats[[3]]-comparison[[3]]}},"\[CapitalDelta]\[LeftAngleBracket]x^2\[RightAngleBracket] = "<>ToStringMod[del[[2]]]<>"\n\[CapitalDelta]\[LeftAngleBracket]x^4\[RightAngleBracket] = "<>ToStringMod[del[[3]]]<>"\n\[LeftAngleBracket]H\[RightAngleBracket] = "<>ToStringMod[del[[1]]]]
reportStatsComparison[f_,H_,comparison_]:=reportStatsComparison[calcStats[f,H],comparison]
Opsubs[Op_,subs_]:=(Op[#]/.subs)&
c[n_]:=If[n==0,1,ToExpression["C"<>ToString[n]]];
dn[n_]:=ToExpression["d"<>ToString[n]];
Hermite[n_,even_:2]:=x^(even n) Exp[-dn[n]x^2];
Cpsi[i_]:=c[i]Hermite[i]
Psis[n_]:=Module[{sum},
sum=Sum[Cpsi[i],{i, 0, n}];
sum
]

MinimizeEresult[Eresult_,pow_]:=NMinimize[Join[{Eresult,Eresult>=0},Table[dn[n]>0,{n,0,pow}],Table[c[n]>=0,{n,1,pow}]],Join[Table[dn[n],{n,0,pow}],Table[c[n],{n,1,pow}]],AccuracyGoal->6]
CalcEnergies[psis_]:=IntInf[#*H[#]]&/@psis;
CalcStats[psis_,energies_,params_]:=Module[{len,enparams,sols,psisnormed},
len=Length[energies];
enparams=energies/.params;
sols=Table[MinimizeEresult[enparams[[i]],i-1],{i,1,len}];
psisnormed=Table[psis[[i]]/.sols[[i]][[2]],{i,1,len}];
Table[{sols[[i]][[2]],Join[{0},exactstats/.params],Join[{sols[[i]][[1]]},calcStatsNoH[psisnormed[[i]]]]},{i,1,len}]
];
GenFigure[psis_,stats_,params_,letter_]:=Module[{psisubs,plots,legends},
psisubs=Table[psis[[i]]/.params/.stats[[i]][[1]],{i,1,Length[psis]}];
plots=Join[{exactnormed/.params},psisubs];
legends=Join[{"Exact"},Table[If[n==1,"Gaussian\n",ToString[n]<>"-Hermite\n"]<>reportStatsComparison[stats[[n]][[3]],stats[[n]][[2]]],{n,1,Length[psis]}]];
Show[{Plot[plots,{x,-5,5},PlotLegends->legends,PlotRange->Full,PlotLabel->"c = "<>ToString[c/.params]<>", \[CapitalGamma] = "<>ToString[\[CapitalGamma]/.params]],letter},ImageSize->All]
];

(* ::Package:: *)

(* ::Title:: *)
(*A Mathematica NURBS (Non-Uniform Rational B-Splines) Package and Its Applications*)


(* ::Text:: *)
(*Author: Liutong Zhou*)


BeginPackage["NURBS`"]
NDegreePowerBasis::usage="NDegreePowerBasis[pts] creates an N Degree Power Basis Curve with control points pts."
Bernstein::usage="Bernstein[i, n, v] returns the ith term of the n-degree Bernstein polynomial basis"
MyBezierCurve::usage="MyBezierCurve[points,u] Returns the parametric form of the BezierCurve definend by the control points"
PlotBezierCurve::usage="PlotBezierCurve[points] plots the BezierCurve defined by the control points"
RationalBezierCurve::usage="RationalBezierCurve[points,u] gives the parametric form of the Rational Bezier Curve definend by the control points. 
The last column of points is the weight of the points"
PlotRationalBezierCurve::usage="PlotRationalBezierCurve[data] plots the RationalBezierCurve[points,u]"
PowerBasisSurface::usage="PowerBasisSurface[points,u,v] gives the Power Basis Surface defined by control points in parametric form S(u,v)"
PlotPowerBasisSurface::usage="PlotPowerBasisSurface plots the PowerBasisSurface defined by control points"
BezierSurface::usage="BezierSurface[points,u,v] gives out the Bezier Surface S(u,v) defined by the control points in Parametric form"
PlotBezierSurface::usage="PlotBezierSurface[points] plot the Bezier Surface S(u,v) defined by the control ponints"
RationalBezierSurface::usage="RationalBezierSurface[points,u,v] gives the parametric form of the Rational Bezier Surface definend by the control points. 
The last column of the points is the weight of the points"
PlotRationalBezierSurface::usage="PlotRationalBezierSurface[points] plots the RationalBezierSurface[points,u,v]"
MyBSplineBasis::usage="MyBsplineBasis[U,i,p,u] returns the ith p degree B-Spline Basis function defined by knotvector U"
MyBSplineCurve::usage="MyBSplineCurve[points,U,p,u] returns the parametric form of the p degree B-Spline Curve defined by knot vector U"
PlotBSplineCurve::usage="PlotBSplineCurve[points,U,p] plots the p degree B Spline curve defined by knot vector U"
DBSplineCurve::usage="DBSplineCurve[points,U,p,u] calculates the first-order derivative of the p degree B-Spline Curve C(u)"
MyBSplineSurface::usage="MyBSplineSurface[points,U,V,p,q,u,v] gives the Parametric function S(u,v) of p,q degree B-Spline Surface defined by control points and knot vectors U, V "
PlotBSplineSurface::usage=" PlotBSplineSurface[points,U,V,p,q] plots the p q degree B-Spline Surface defined by cotrol points and knot vectors U, V"
Begin["Private`"]


(* ::Chapter::Closed:: *)
(*Chapter 1 Curve and Surface Basics*)


(* ::Section:: *)
(*1.2 Power Basis Form of a Curve*)


NDegreePowerBasis[points:{{_, _} ..} | {{_, _, _} ..}]:=Module[{l=Length[points],dimension=Length@points[[1]],parametricform,transpts,u},
transpts=Transpose[points];
parametricform=transpts.Table[u^n,{n,1,l}];
If [dimension==3,ParametricPlot3D[parametricform,{u,0,1},PlotRange->MinMax/@transpts,AspectRatio->1],ParametricPlot[parametricform,{u,0,1},PlotRange->MinMax/@transpts,AspectRatio->1]]
]


(* ::Section:: *)
(*1.3 Bezier Curves*)


Bernstein[i_, n_, v_] := Binomial[n, i] v^i*(1 - v)^(n - i) ;

MyBezierCurve[pt : {{_, _} ..} | {{_, _, _} ..}, v_] :=Module[{n = Length[pt] - 1},
  Simplify@Sum[Bernstein[i, n, v]*pt[[i+1]], {i, 0, n}]
];

PlotBezierCurve[pt : {{_, _} ..} | {{_, _, _} ..}] := 
    Module[{BC,v,CurveGraph,CurvePolygon,CurvePoints},
      BC = MyBezierCurve[pt, v];
      If[Length[pt[[1]]] == 2, Module[{},
                 CurveGraph = ParametricPlot[BC, {v, 0, 1}, PlotRange -> All];
                 CurvePolygon = ListPlot[pt, Joined -> True];
                 CurvePoints =  ListPlot[pt, Joined -> False, PlotStyle -> PointSize[0.02]]
               ], 
           Module[{},
                     CurveGraph = ParametricPlot3D[BC, {v,0, 1},PlotRange -> All];
                     CurvePolygon = Graphics3D[Line[pt]];
                     CurvePoints = Graphics3D[Map[{PointSize[0.02], Point[#]} &, pt]]
          ]];
      Show[CurveGraph, CurvePolygon, CurvePoints]];


(* ::Section:: *)
(*1.4 Rational Bezier Curve*)


RationalBezierCurve[pts_,u_]:=Module[{n=Length[pts]-1,i},
(Table[Bernstein[i,n,u],{i,0,n}].(pts[[All,1;;-2]]*pts[[All,-1]]))/
(Table[Bernstein[i,n,u],{i,0,n}].pts[[All,-1]])//Simplify];

PlotRationalBezierCurve[pt : {{_, _, _} ..} | {{_, _, _, _} ..}] := 
    Module[{RBC,v,CurveGraph,CurvePolygon,CurvePoints},    
    RBC=  RationalBezierCurve[pt, v];
    If[Length[pt[[1]]] == 3, Module[{},   
 CurveGraph =ParametricPlot[RBC, {v, 0, 1}];
CurvePolygon =ListPlot[pt[[All,1;;2]],Joined->True,PlotStyle->Black];
CurvePoints = ListPlot[pt[[All,1;;2]],PlotStyle->{Black, PointSize[0.02]},PlotRange->All];] 
  , Module[{},         
CurveGraph =   ParametricPlot3D[   RBC, {v, 0,  1}];     
CurvePolygon =         Graphics3D[Line[pt[[All,1;;-2]]] ];
CurvePoints = Graphics3D[Point[pt[[All,1;;-2]]],PlotRange->All];]
]
CurveGraph
Show[  CurvePoints,CurvePolygon ,CurveGraph]];


(* ::Section:: *)
(*1.5 Tensor Product Surfaces*)


(* ::Text:: *)
(*PowerBasisSurface*)


PowerBasisSurface[points:{{{_, _, _} ..}..},u_,v_]:=Module[{d=Dimensions[points,2]-1,f,g,n,m},
n=d[[1]];   m=d[[2]];
{f, g} =   MapThread[Table[#2^i, {i, 0, #1}] &, {{n, m}, {u, v}}];
f.points*g//Total//Simplify
]

PlotPowerBasisSurface[points_]:=Module[{u,v},
Show[Graphics3D[{PointSize[0.015],Point[#]}&@Flatten[points,1]],ParametricPlot3D[PowerBasisSurface[points,u,v],{u,0,1},{v,0,1}]]
]


(* ::Text:: *)
(*Bezier Surface*)


BezierSurface[points : {{{_, _, _} ..} ..}, u_, v_] :=     Module[{d=Dimensions[points,2]-1,n,m, f,g},
n=d[[1]];   m=d[[2]];
      {f, g} =        MapThread[  Table[Bernstein[i, #1, #2], {i, 0, #1}] &, {{n, m}, {u, v}}];
f.points*g//Total//Simplify];

PlotBezierSurface[pts : {{{_, _, _} ..} ..}] := 
    Module[{BS,u,v,SurfaceGraph,SurfacePoints},
      BS = BezierSurface[pts, u,v];
      SurfaceGraph =         ParametricPlot3D[BS, {u, 0, 1}, {v, 0, 1}];
      SurfacePoints =         Graphics3D[{PointSize[0.015], Point[#]} &@ Flatten[pts, 1]];
      Show[SurfacePoints,SurfaceGraph]];


(* ::Text:: *)
(*Rational Bezier Surface*)


RationalBezierSurface[points : {{{_, _, _, _} ..} ..}, u_, v_] :=     Module[{d=Dimensions[points,2]-1, n,m,f, g,sumweight,R},
n=d[[1]];m=d[[2]];
      {f, g} =MapThread[Table[Bernstein[i, #1, #2], {i, 0, #1}] &, {{n, m}, {u, v}}];
sumweight=f.points[[;;,;;,-1]].g;
R=Outer[Times,f,g]*points[[;;,;;,-1]]/sumweight;
Total[R*points[[;;,;;,1;;-2]],2]//Simplify];

PlotRationalBezierSurface[pts : {{{_, _, _, _} ..} ..}] :=     Module[{u,v,RBS,SurfaceGraph,SurfacePoints},
      RBS = RationalBezierSurface[pts, u,v];
      SurfaceGraph = ParametricPlot3D[RBS, {u, 0, 1}, {v, 0, 1}];
      SurfacePoints =        Graphics3D[{PointSize[0.015], Point[#]} &@ Flatten[pts[[;;,;;,1;;-2]],1]];
      Show[  SurfacePoints,SurfaceGraph]];


(* ::Chapter:: *)
(*Chapter 2  B-Spline Basis Function*)


(* ::Section:: *)
(*2.2 B-Spline Basis Function*)


MyBSplineBasis[U_List,i_,0,u_]:=If[U[[i+1]]<= u<U[[i+2]],1,0];
MyBSplineBasis[U_List,i_,p_,u_]:=If[MyBSplineBasis[U,i,p-1,u]==0,0,(u-U[[i+1]])/(U[[i+1+p]]-U[[i+1]])*MyBSplineBasis[U,i,p-1,u]]+
If[MyBSplineBasis[U,i+1,p-1,u]==0,0,(U[[i+1+p+1]]-u)/(U[[i+1+p+1]]-U[[i+1+1]])* MyBSplineBasis[U,i+1,p-1,u]]//Simplify;


(* ::Chapter:: *)
(*Chapter 3 B-Spline Curves and Surfaces*)


(* ::Section:: *)
(*3.2 B-Spline Curve*)


MyBSplineCurve[points : {{_, _} ..} | {{_, _, _} ..}, U_List, p_Integer,u_] := Module[{n=Length[points]-1},
Sum[points[[i+1]]*MyBSplineBasis[U,i,p,u],{i,0,n}]//Simplify
    ];


PlotBSplineCurve[points : {{_, _} ..} | {{_, _, _} ..}, U_, p_] :=     Module[{BSC,u,CurveGraph,ControlPolygon},      
BSC = MyBSplineCurve[points, U, p,u];
      If[Length[points[[1]]] == 2,      
      CurveGraph =  ParametricPlot[BSC,{u,0,1-0.0001}];
      ControlPolygon=Graphics[{{Black,PointSize[0.018],Point[points]},{Green,Thin,Line[points]}},Axes->True]       ,
      CurveGraph=ParametricPlot3D[BSC,{u,0,1-0.0001}];
      ControlPolygon=Graphics3D[{{Black,PointSize[0.018],Point[points]},{Green,Thin,Line[points]}},Axes->True]]; 
Show[ControlPolygon,CurveGraph]
];


(* ::Section:: *)
(*3.3 Derivative of B-Spline Curve*)


DBSplineCurve[points_,U_,p_,u_]:=Sum[
MyBSplineBasis[U,i,p-1,u]*p*(points[[i+1+1]]-points[[i+1]])/(U[[i+1+p+1]]-U[[i+1+1]]),
{i,0,Length[points]-1-1}]//Simplify;


(* ::Section:: *)
(*3.4 B-Spline Surface*)


MyBSplineSurface[points:{{{_,_,_}..}..},U_,V_,p_,q_,u_,v_]:=Module[{d=Dimensions[points],n,m},
n=d[[1]]-1;m=d[[2]]-1;
Sum[MyBSplineBasis[U,i,p,u]*MyBSplineBasis[V,j,q,v]*points[[i+1,j+1]],{i,0,n},{j,0,m}]//Simplify
]

PlotBSplineSurface[points:{{{_,_,_}..}..},U_,V_,p_,q_]:=Show[Graphics3D[{Black,PointSize[0.02],Point[Flatten[points,1]]}],
ParametricPlot3D[BSplineFunction[points,SplineDegree->{p,q}][u,v],{u,0,1},{v,0,1}]
];


End[]
EndPackage[]
Print["Package NURBS.wl By Liutong Zhou, Version 1.0\nCreated on Nov. 25 2016\nLast Updated on Dec 11 2016\n
User-accessible functions:
NDegreePowerBasis,Bernstein, MyBezierCurve, PlotBezierCurve, RationalBezierCurve, PlotRationalBezierCurve,BezierSurface,PlotBezierSurface,
RationalBezierSurface,PlotRationalBezierSurface, MyBSplineBasis, MyBSplineCurve, PlotBSplineCurve, DBSplineCurve, MyBSplineSurface, PlotBSplineSurface\n
Type ?functionname to see usage information
"]

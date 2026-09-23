(* n=3 tree EYM / YM collinear kernel, exact *)
ClearAll["Global`*"];
eta = DiagonalMatrix[{1, -1, -1, -1}];
lower[v_] := eta . v;
Dot4[a_, b_] := Dot[lower[a], b];
sig0 = {{1, 0}, {0, 1}};
sig1 = {{0, 1}, {1, 0}};
sig2 = {{0, -I}, {I, 0}};
sig3 = {{1, 0}, {0, -1}};
sigma = {sig0, sig1, sig2, sig3};
sigbar = {sig0, -sig1, -sig2, -sig3};
Ang[i_, j_] := -(lam[i][[1]]*lam[j][[2]] - lam[i][[2]]*lam[j][[1]]);
Sq[i_, j_] := lamT[i][[1]]*lamT[j][[2]] - lamT[i][[2]]*lamT[j][[1]];
Momentum[i_] := Module[{M},
  M = Outer[Times, lam[i], lamT[i]];
  Table[(1/2) Tr[sigbar[[mu]] . M], {mu, 4}]
];
AngSqVec[i_, j_] := Module[{M},
  M = Outer[Times, lam[i], lamT[j]];
  Table[Tr[sigbar[[mu]] . M], {mu, 4}]
];
SqAngVec[i_, j_] := Module[{M},
  M = Outer[Times, lam[j], lamT[i]];
  Table[Tr[sigbar[[mu]] . M], {mu, 4}]
];
ePlus[k_, r_] := AngSqVec[r, k]/(Sqrt[2] Ang[r, k]);
eMinus[k_, r_] := -SqAngVec[r, k]/(Sqrt[2] Sq[r, k]);
Kin3[e1_, e2_, e3_, p_, q_, r_] :=
  Dot4[e1, e2]*Dot4[p - q, e3] + Dot4[e2, e3]*Dot4[q - r, e1] + Dot4[e3, e1]*Dot4[r - p, e2];
eps3 = LeviCivitaTensor[3];
f2[a_, b_, c_] := eps3[[a, b, c]];

(* abelian graviton-two-gluon vertex, i*L cross linear in h *)
Hgg[p1_, e1_, p2_, e2_, hlow_] := Module[
  {F1, F2, htrace, hup, Lonly, crossZ},
  F1[mu_, nu_] := -I*(lower[p1][[mu]]*lower[e1][[nu]] - lower[p1][[nu]]*lower[e1][[mu]]);
  F2[mu_, nu_] := -I*(lower[p2][[mu]]*lower[e2][[nu]] - lower[p2][[nu]]*lower[e2][[mu]]);
  htrace = Tr[eta . hlow];
  hup = eta . hlow . eta;
  Lonly[z_, u1_, u2_] := Module[{ginv, sqrtg, Fa, Fb, Lc, mu, nu, aa, bb},
    ginv = eta - z*hup;
    sqrtg = 1 + z*htrace/2;
    Lc = 0;
    Do[
      Fa = If[u1, F1[mu, nu], 0] + If[u2, F2[mu, nu], 0];
      Fb = If[u1, F1[aa, bb], 0] + If[u2, F2[aa, bb], 0];
      Lc += -1/4*sqrtg*ginv[[mu, aa]]*ginv[[nu, bb]]*Fa*Fb;
      , {mu, 4}, {nu, 4}, {aa, 4}, {bb, 4}];
    Lc
  ];
  crossZ = (Lonly[1, True, True] - Lonly[0, True, True]) -
    (Lonly[1, True, False] - Lonly[0, True, False]) -
    (Lonly[1, False, True] - Lonly[0, False, True]);
  I*crossZ
];

(* full contact+exchange for colors, Lagrangian normalization *)
EYM[cols_, pols_, moms_, eG_, pG_] := Module[
  {hmn, contact, exch, jj, ii, kk, q, a, vh, vg, basis, s, hlow},
  hlow = Outer[Times, lower[eG], lower[eG]];
  hmn = hlow;
  basis = IdentityMatrix[4];
  (* contact: use Hgg-style by calling a dedicated cubic routine below via global ContactEYM *)
  contact = ContactEYM[moms, pols, cols, hlow];
  exch = 0;
  Do[
    ii = Mod[jj, 3] + 1;
    kk = Mod[jj + 1, 3] + 1;
    q = pG + moms[[jj]];
    s = Dot4[q, q];
    vh = Table[Hgg[moms[[jj]], pols[[jj]], -q, basis[[a]], hmn], {a, 4}];
    vg = Table[
      f2[cols[[ii]], cols[[kk]], cols[[jj]]]*
       Kin3[pols[[ii]], pols[[kk]], basis[[a]], moms[[ii]], moms[[kk]], q],
      {a, 4}];
    exch += (-I)/s*Sum[vh[[a]]*eta[[a, a]]*vg[[a]], {a, 4}];
    , {jj, 3}];
  contact + exch
];

ContactEYM[gMoms_, gPols_, gCols_, hlow_] := Module[
  {zg, zh, F, L, mu, nu, a, aa, bb, expr, hup, htrace, ginv, sqrtg, b, c, k, Ab, Ac, n = 3},
  htrace = Tr[eta . hlow];
  hup = eta . hlow . eta;
  ginv = eta - zh*hup;
  sqrtg = 1 + zh*htrace/2;
  F = Table[0, {a, 3}, {mu, 4}, {nu, 4}];
  Do[
    F[[gCols[[k]], mu, nu]] += zg[k]*(-I)*(lower[gMoms[[k]]][[mu]]*lower[gPols[[k]]][[nu]] -
         lower[gMoms[[k]]][[nu]]*lower[gPols[[k]]][[mu]]);
    , {k, n}, {mu, 4}, {nu, 4}];
  Do[
    Ab = Sum[zg[k]*lower[gPols[[k]]][[mu]]*Boole[gCols[[k]] == b], {k, n}];
    Ac = Sum[zg[k]*lower[gPols[[k]]][[nu]]*Boole[gCols[[k]] == c], {k, n}];
    F[[a, mu, nu]] += f2[a, b, c]*Ab*Ac;
    , {a, 3}, {b, 3}, {c, 3}, {mu, 4}, {nu, 4}];
  L = Sum[-1/4*sqrtg*ginv[[mu, aa]]*ginv[[nu, bb]]*F[[a, mu, nu]]*F[[a, aa, bb]],
    {a, 3}, {mu, 4}, {nu, 4}, {aa, 4}, {bb, 4}];
  expr = Expand[I*L];
  Coefficient[expr, zh*zg[1]*zg[2]*zg[3]]
];

SetSpinors[Ls_, LTs_] := Do[lam[i] = Ls[[i]]; lamT[i] = LTs[[i]];, {i, Length[Ls]}];

Solve4[L1_, LT1_, L2_, LT2_, L3_, L4_] := Module[{sol},
  sol = Solve[
    Flatten[
      Outer[Times, L1, LT1] + Outer[Times, L2, LT2] +
        Outer[Times, L3, {a3, b3}] + Outer[Times, L4, {a4, b4}]
    ] == {0, 0, 0, 0},
    {a3, b3, a4, b4}
  ][[1]];
  {L1, LT1, L2, LT2, L3, {a3, b3} /. sol, L4, {a4, b4} /. sol}
];

AmpStandard[Ls_, LTs_, ref_] := Module[{moms, pols, mf},
  SetSpinors[Join[Partition[Riffle[Ls, LTs], 2][[All, 1]] (* dummy *), {}], {}];
  (* set explicitly *)
  Do[lam[i] = Ls[[i]]; lamT[i] = LTs[[i]];, {i, 4}];
  lam[5] = ref[[1]]; lamT[5] = ref[[2]];
  moms = Table[Momentum[i], {i, 4}];
  pols = {eMinus[1, 5], eMinus[2, 5], ePlus[3, 5], ePlus[4, 5]};
  mf = EYM[{1, 2, 3}, pols[[1 ;; 3]], moms[[1 ;; 3]], pols[[4]], moms[[4]]];
  {Simplify[mf/Sqrt[2]], Simplify[Total[moms]]}
];

(* reference point *)
pt0 = Solve4[{1, 0}, {1, 2}, {0, 1}, {3, -1}, {1, 1}, {2, -1}];
Ls0 = {pt0[[1]], pt0[[3]], pt0[[5]], pt0[[7]]};
LTs0 = {pt0[[2]], pt0[[4]], pt0[[6]], pt0[[8]]};
res0 = AmpStandard[Ls0, LTs0, {{1, 3}, {2, 5}}];
Print["POINT0 ", res0[[1]], " cons ", res0[[2]]];

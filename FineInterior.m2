-- -*- coding: utf-8 -*-
newPackage(
    "FineInterior",
    Version => "0.1",
    Date => "June 8, 2023",
    Authors => {
	{Name => "Simen Moe", Email => "swm418@ic.ac.uk"}},
    Headline => "Compute lattice width and Fine interior of lattice polytopes.",
    Keywords => {"Documentation"},
    DebuggingMode => false,
    PackageImports=>{"NormalToricVarieties", "MonomialAlgebras",  "Polyhedra", "ReesAlgebra", "MinimalPrimes", "PrimaryDecomposition"}
    )

export {"latticeWidthCoord", "latticeWidth", "ord", "fineInterior"}

read 

ord = method(TypicalValue => QQ)
ord (Polyhedron, Matrix) := (P,n) -> (
    Ppts := latticePoints P;
    dists := for xx in Ppts list (((transpose(xx)*n)_0)_0);
    min(dists))

fineInterior = method(TypicalValue => QQ)
fineInterior Matrix := M -> (
    P := convexHull(M);
    m := dim P;
    vertP := vertices P;
    Nlist := for d from 1 to m list (
    	L := faces(d,P);
	dFaces := apply(L, f -> vertP_(f#0));
	dbasis := for dface in dFaces list (
    	    Q := convexHull(matrix dface);	 
    	    nCone := normalCone(P,Q);
	    hBas := hilbertBasis nCone;
	    hBas);
     	 dbasis);
    Nlist = flatten flatten Nlist;
    Nlist = unique Nlist;
    Nlist = Nlist;
    v2 := for w2 in Nlist list (-ord(P,w2)-1);
    Nlist = (-1)*Nlist;
    v2 = matrix {v2};	
    v2 = transpose v2;
    M2 :=  matrix {Nlist};
    M2 = transpose M2;
    P2 := polyhedronFromHData(M2,v2);
    P2)

fineInterior Polyhedron := P -> (
    M := vertices P;
    fineInterior(M))

latticeWidth = method()
latticeWidth (Polyhedron, ZZ) := (P,maxWidth) -> (
    if not isFullDimensional(P) then (error("Error: Not full dimensional"));
    minWidth := maxWidth;
    m := dim P;
    vertP := vertices P;
    Nlist := for d from 1 to m list (
    	L := faces(d,P);
	    dFaces := apply(L, f -> vertP_(f#0));
        dbasis := for dface in dFaces list (
                    Q := convexHull(matrix dface);	 
                    nCone := normalCone(P,Q);
                    hBas := hilbertBasis nCone;
                    hBas);
     	dbasis);
    Nlist = flatten flatten Nlist;
    Nlist = unique Nlist;
    Nlist = (-1)*Nlist;
    Nlist = matrix{Nlist};
    MM := transpose(Nlist)*vertP;
    minVal := maxWidth;
    numdiffs := for i from 0 to rank target MM-1 list (if minVal==1 then break else (ww := flatten entries MM^{i}; minValTemp := max(ww)-min(ww);if minValTemp < minVal then minVal=minValTemp));
    if (minVal==maxWidth) then print("Warning: Maxwidth equals minwidth. Result may be inaccurate.");
    minVal
)

latticeWidth (Matrix, ZZ) := (M,maxWidth) -> (
    P := convexHull M;
    latticeWidth(P,maxWidth)
)

latticeWidth (Matrix) := M -> (
    latticeWidth(M,5)
)

latticeWidth (Polyhedron) := P -> (
    latticeWidth(P,5)
)

latticeWidthCoord = method()
latticeWidthCoord (Polyhedron, ZZ) := (P,maxWidth) -> (
    if not isFullDimensional(P) then (error("Warning: Not full dimensional"));
    finalPlane := 0;
    minWidth := maxWidth;
    m := dim P;
    vertP := vertices P;
    Nlist := for d from 1 to m list (
    	L := faces(d,P);
	dFaces := apply(L, f -> vertP_(f#0));
	dbasis := for dface in dFaces list (
    	    Q := convexHull(matrix dface);	 
    	    nCone := normalCone(P,Q);
	    hBas := hilbertBasis nCone;
	    hBas);
     	 dbasis);
    Nlist = flatten flatten Nlist;
    Nlist = unique Nlist;
    Nlist = (-1)*Nlist;
    Nlist = matrix{Nlist};
    MM := transpose(Nlist)*vertP;
    minVal := maxWidth;
    numdiffs := for i from 0 to rank target MM-1 list (ww := flatten entries MM^{i}; minValTemp := max(ww)-min(ww);if minValTemp < minVal then (minVal=minValTemp;finalPlane=Nlist_{i}));
    if (minVal==maxWidth) then print("Warning: Maxwidth equals minwidth. Results may be inaccurate.");
    (minVal,finalPlane)
)

latticeWidthCoordAll = method()
latticeWidthCoordAll (Polyhedron, ZZ) := (P,maxWidth) -> (
    if not isFullDimensional(P) then (error("Warning: Not full dimensional"));
    finalPlanes := {};
    minWidth := maxWidth;
    m := dim P;
    vertP := vertices P;
    Nlist := for d from 1 to m list (
    	L := faces(d,P);
	dFaces := apply(L, f -> vertP_(f#0));
	dbasis := for dface in dFaces list (
    	    Q := convexHull(matrix dface);	 
    	    nCone := normalCone(P,Q);
	    hBas := hilbertBasis nCone;
	    hBas);
     	 dbasis);
    Nlist = flatten flatten Nlist;
    Nlist = unique Nlist;
    Nlist = (-1)*Nlist;
    Nlist = matrix{Nlist};
    MM := transpose(Nlist)*vertP;
    minVal := maxWidth;
    numdiffs := for i from 0 to rank target MM-1 list (ww := flatten entries MM^{i}; minValTemp := max(ww)-min(ww);if minValTemp <= minVal then (minVal=minValTemp;finalPlanes=append(finalPlanes,[minVal,Nlist_{i}])));
    if (minVal==maxWidth) then print("Warning: Maxwidth equals minwidth. Results may be inaccurate.");
    finalPlanes2 := for pair in finalPlanes list (
        if pair_0==minVal then pair else continue
    );
    (minVal,finalPlanes2)
)

latticeWidthCoord (Matrix, ZZ) := (M,maxWidth) -> (
    P := convexHull M;
    latticeWidthCoord(P,maxWidth)
)

latticeWidthCoord (Matrix) := M -> (
    latticeWidthCoord(M,5)
)

latticeWidthCoord (Polyhedron) := P -> (
    latticeWidthCoord(P,5)
)

latticeWidthCoordAll (Matrix, ZZ) := (M,maxWidth) -> (
    P := convexHull M;
    latticeWidthCoordAll(P,maxWidth)
)

latticeWidthCoordAll (Matrix) := M -> (
    latticeWidthCoordAll(M,5)
)

latticeWidthCoordAll (Polyhedron) := P -> (
    latticeWidthCoordAll(P,5)
)


--------------------
-- DOCUMENTATION 
--------------------

beginDocumentation()

document {
     Key => {latticeWidthCoord, (latticeWidthCoord,Polyhedron)},
     Headline => "Computes the lattice width of a lattice polytope and gives a hyperplane obtaining the width",
     Usage => "(lw,p) = latticeWidthCoord M",
     Inputs => {
	  "P" => Polyhedron
	  },
     Outputs => {
	  "lw" => ZZ,
      "p" => Matrix
	  },

     PARA{}, TT "latticeWidthCoord", "Computes the lattice width of the polytope P. Returns the lattice width and a hyperplane to which this width is obtained."
     }

document {
     Key => {latticeWidth, (latticeWidth,Polyhedron)},
     Headline => "Compute the lattice width of lattice polytope",
     Usage => "lw = latticeWidth P",
     Inputs => {
	  "P" => Polyhedron
	  },
     Outputs => {
	  "lw" => ZZ
	  },

     PARA{}, TT "latticeWidth", "Compute the lattice width of a lattice polytope."
     }


document {
     Key => {fineInterior, (fineInterior,Polyhedron)},
     Headline => "Computes the Fine interior.",
     Usage => " PFI = fineInterior P",
     Inputs => {
	  "P" => Polyhedron
	  },
     Outputs => {
	  "PFI" => Polyhedron
	  },

     PARA{}, TT "fineInterior", "Computes the Fine interior of a lattice polytope."
     }
end


restart
path = append(path, "/Users/simenwm/Dropbox/Projects/Macaulay2/Macaulay2_Packages/");
stack path
uninstallPackage "FineInterior"
installPackage "FineInterior"

needsPackage "Polyhedra"
M = transpose matrix{{0,2,2},{1,3,0},{2,4,3},{3,0,1}}
P = convexHull M
latticeWidth P
latticeWidthCoord P
PFI = fineInterior P
dim PFI
vertices PFI
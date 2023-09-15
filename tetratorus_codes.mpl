restart:
libname:=libname, "/home/data/archiv/daniel/maple/lib10": 
libname := libname, "/users/kirekod/maple/lib":
with(Involutive):with(LinearAlgebra):
with(SimplicialSurfaceEmbeddings):
with(combinat):



AddTetra:=proc(surf,i,vert,f)
local g,c,v,v1,v2,v3,n,sol,vv;
v1:=CoordinateMatrix(surf,i,"vertices"=[f[2]])-CoordinateMatrix(surf,i,"vertices"=[f[1]]);
v:=CoordinateMatrix(surf,i,"vertices"=[f[1]]);
v2:=CoordinateMatrix(surf,i,"vertices"=[f[3]])-CoordinateMatrix(surf,i,"vertices"=[f[1]]);
n:=CrossProduct(v1,v2);v3:=CoordinateMatrix(surf,i,"vertices"=[vert]);
vv:=[v3[1][1]+t*n[1],v3[2][1]+t*n[2],v3[3][1]+t*n[3]];
sol:=solve(n[1]*vv[1]+n[2]*vv[2]+n[3]*vv[3]=n[1]*v[1][1]+n[2]*v[2][1]+n[3]*v[3][1],{t});
return subs(sol,[v3[1][1]+2*t*n[1],v3[2][1]+2*t*n[2],v3[3][1]+2*t*n[3]]);
end proc:



AddT:=proc(surface,face,V)
local g,newV,vertices,c,vof,copyVOF,newVertex,s,l;
RemoveFace(surface,convert(face,list));
newV:=max(Vertices(surface))+1;vertices:=[op(Vertices(surface)),newV];
c:=CoordinateMatrix(surface,1,"listlist"=true);vof:=Faces(surface);
vof:={seq({vof[l][1],vof[l][2],vof[l][3]},l=1..nops(vof))}; 
vof:=(vof minus {face}) union {{face[1],face[2],newV},{face[2],face[3],newV},{face[1],face[3],newV}};
copyVOF:=[seq([vof[l][1],vof[l][2],vof[l][3]],l=1..nops(vof))];
newVertex:=AddTetra(surface,1,V,[face[1],face[2],face[3]]);
c:=[op(c),newVertex];
s:=NewSurface();DefineEmbedding(s,c,"vertices"=vertices,"faces"=copyVOF);
return s;
end proc:



# given our "symbol" cacsym for a cactus ( cacsym=[[face,vertex],...]) we want to construct a cacti by iteratively reflecting the vertex 'vertex' through the plane spanned by the face 'face' 
ConstructCactus:=proc(symbol)
local g,c,x,h,surf;x:=a;h:=b;
c:=[[1/2, 0, 0], [-1/2, 0, 0], [(x^2 - 1)/(2*(x^2 + 1)), x/(x^2 + 1), h], [-(x^2 - 1)/(2*(x^2 + 1)), -x/(x^2 + 1), h]];
surf:=NewSurface(): DefineEmbedding(surf,c,"faces"=[[1,2,3],[1,2,4],[1,3,4],[2,3,4]],"vertices"=[1,2,3,4]);
for g in symbol do 
  surf:=AddT(surf,g[1],g[2]);
end do;
return surf;
end proc: 



DegreesOfVertices:=proc(surface)
local g,i,f,degrees,count; degrees:=[];
for i in Vertices(surface) do 
  count:=0;
  for f in Faces(surface) do
    if member(i,f) then 
      count:=count+1;
    end if;
  end do;
  degrees:=[op(degrees), count];
end do;
return degrees;
end proc:



NeighbourVertices:=proc(surface,v)
local g,res,edges,e;
edges:=Edges(surface);res:=[];
for e in Edges do
  if e[1]=v then 
    res:=[op(res),e[2]] 
  end if;
  if e[2]=v then 
    res:=[op(res),e[1]] 
  end if;
od; 
return res;
end proc:



FacesOfVertex:=proc(surface,v)
local g,res,faces,f;
res:=[];faces:=Faces(surface);
for f in faces do
  if f[1]=v or f[2]=v or f[3]=v then 
    res:=[op(res),f]; 
  end if;
end do; 
return res;
end proc:



help_norm:=proc(v)
return sqrt(v[1]^2+v[2]^2+v[3]^2);
end proc:



help_check := proc(surface, face1, face2) 
local c, v11, v12, v13, v21, v22, v23, d1, d2, d3; 
c := simplify(CoordinateMatrix(surface, 1, "listlist" = true)); 
v11 := c[face1[1]]; v12 := c[face1[2]]; v13 := c[face1[3]]; 
v21 := c[face2[1]]; v22 := c[face2[2]]; v23 := c[face2[3]]; 
d1 := Determinant(<<v12 - v11> | <v13 - v11> | <v21 - v11>>); 
d2 := Determinant(<<v12 - v11> | <v13 - v11> | <v22 - v11>>); 
d3 := Determinant(<<v12 - v11> | <v13 - v11> | <v23 - v11>>); 
return [simplify(d1), simplify(d2), simplify(d3)]; 
end proc:



VerticesOfDegreeThree:=proc(surface)
local g,res,v;
res:=[];
for v in Vertices(surface) do 
  if DegreesOfVertices(surface)[v]=3 then 
    res:=[op(res),v]; 
  end if;
end do;
return res;
end proc:



help_norm:=proc(v)
return sqrt(v[1]^2+v[2]^2+v[3]^2);
end proc:

IsVertexfaithful:=proc(coordinates)
local g,temp,i,coor;
coor:=coordinates;
for i from 1 to nops(coor) do 
  temp:=select(v->evalf(help_norm(coor[i]-v))<=10.0^(-10),coor); 
  if nops(temp)>1 then 
    return false;  
  end if;
end do;
return true;
end proc:



IsVertexfaithful_OnlyFaces:=proc(coordinates,face1,face2) # only face1 and face2 should have coordinates of vertices in common and other coordinates should be vertex faithful. 
local g,temp,i,coor,remaining_vertices,temp1;
coor:=coordinates;
temp:=[]; temp1:=[];
remaining_vertices:={op([$1..nops(coor)])} minus {op(face1),op(face2)}; 
for i in remaining_vertices do 
  temp:=map(v->evalf(help_norm(coor[i]-v)),coor);    #before 10.0^(-10)
  temp1:=select(t->10.0^(-8)>=t,temp);
  if nops(temp1)>1 then
    return false;  
  end if;
end do;
temp:=[]; 
temp1:=[];
for i in face1 do 
  temp:=map(j->evalf(help_norm(coor[i]-coor[j])),face2); 
  # temp1:=select(t->10.0^(-8)>=t,temp);
  if nops(temp1)<1 then 
    return false;  
  end if;
end do;
return true;
end proc:



HasSelfIntersections1:=proc(s,coor)
local g,c,f,e,v,v1,v2,w,w1,mat,sol;
c:=coor;
for f in Faces(s) do
  for e in select(e-> not e[1] in f and not e[2] in f,Edges(s)) do  
    v:=c[f[1]];v1:=c[f[2]]-c[f[1]]; v2:=c[f[3]]-c[f[1]];
    w:=c[e[1]];w1:=c[e[2]]-c[e[1]];
    mat:=<<v1>|<v2>|<-w1>>; 
    if evalf(Determinant(mat)^2)>= 10.^(-10) then 
      sol:= MatrixInverse(mat).<w-v>;
      if evalf(sol[3])<= 1. and evalf(sol[3])>=0. and evalf(sol[2])<= 1. and evalf(sol[2])>=0. and evalf(sol[1])<= 1. and evalf(sol[1])>=0. and evalf(sol[1]+sol[2])<=1. then 
        return true; 
      end if;
    end if;
  end do;
end do;
return false;
end proc:



HasSelfIntersections2:=proc(s,coor,face1,face2,id)
local g,c,f,e,v,v1,v2,w,w1,mat,sol,temp;
c:=coor;
if id=1 then
  for f in Faces(s) do
    for e in select(e-> not e[1] in f and not e[2] in f,Edges(s)) do  
      v:=c[f[1]];v1:=c[f[2]]-c[f[1]]; v2:=c[f[3]]-c[f[1]];
      w:=c[e[1]];w1:=c[e[2]]-c[e[1]];
      mat:=<<v1>|<v2>|<-w1>>; 
      if evalf(Determinant(mat)^2)>= 10.^(-10) then 
        sol:= MatrixInverse(mat).<w-v>;
        if evalf(sol[3])<= 1. and evalf(sol[3])>=0. and evalf(sol[2])<= 1. and evalf(sol[2])>=0. and evalf(sol[1])<= 1. and evalf(sol[1])>=0. and evalf(sol[1]+sol[2])<=1. then 
          return true; 
        end if;
      end if;
    end do;
  end do;
else                       # without mirror symmetry
  for f in Faces(s) do 
    if f=face1 or f=face2 then
      temp:=select(e->(not e[1] in face1 and not e[2] in face1) or (not e[1] in face2 and not e[2] in face2),Edges(s));
    else
      temp:=select(e-> not e[1] in f and not e[2] in f,Edges(s));
    end if;
    for e in temp do 
      v:=c[f[1]];v1:=c[f[2]]-c[f[1]]; v2:=c[f[3]]-c[f[1]];
      w:=c[e[1]]; w1:=c[e[2]]-c[e[1]];
      mat:=<<v1>|<v2>|<-w1>>; 
      if evalf(Determinant(mat)^2)>= 10.^(-10) then 
        sol:= MatrixInverse(mat).<w-v>;
        if evalf(sol[3])<= 1. and evalf(sol[3])>=0. and evalf(sol[2])<= 1. and evalf(sol[2])>=0. and evalf(sol[1])<= 1. and evalf(sol[1])>=0. and evalf(sol[1]+sol[2])<=1. then 
          return true; 
        end if;
      end if;
    end do;
  end do;
end if;
return false;
end proc:



check_IdentifiedFaces:=proc(coor,face1,face2)
local g,i,c;
c:=coor;
for i in [1,2,3] do
  if nops(selected(cc-> evalf(help_norm(cc-c[face1[i]]))<=10.^(-5),map(j->c[face2[j]],[1,2,3])))<>1 then 
    return false;
  fi;
end do;
return true;
end proc:



HasVertexInPlane:=proc(surface,face1,face2,coor)
local i,j,v1,v2,v,n,verts,edges,kk,k,eqn_of_plane,ans,w,M,eq,count,c;
c:=coor;
count:=0;
v:=c[face1[1]]; v1:=c[face1[2]]-v; v2:=c[face1[3]]-v;
n:=CrossProduct(v1,v2); 
w:=x*n[1] + y*n[2] + z*n[3]+k; 
eqn_of_plane:=subs(k=-(v[1]*n[1]+v[2]*n[2]+v[3]*n[3]),w); 
verts:={op(Vertices(surface))} minus {op(face1),op(face2)}; 
for i from 1 to nops(verts) do 
  eq:=subs([x=c[verts[i]][1],y=c[verts[i]][2],z=c[verts[i]][3]],eqn_of_plane); 
  if evalf(eq)<10^(-5) and -10^(-5)<evalf(eq)  then 
    count:=count+1; 
  end if; 
end do; 
if count>=1 then 
  return true;
else 
  return false; 
end if;
end proc: 


AreAllVerticesOnOneSideOfPlane := proc(surface, face1, face2, coor) 
local i, v1, v2, v, n, verts, edges, ans, w, M, eq, c, temp, count1, count2, j; 
c := evala(convert(coor, radical)); 
count1 := 0; count2 := 0; temp := []; 
v := c[face1[1]]; 
v1 := c[face1[2]] - v; v2 := c[face1[3]] - v; 
n := CrossProduct(v1, v2); 
verts := {op(Vertices(surface))} minus {op(face1), op(face2)}; 
for j to nops(verts) do 
  temp := [op(temp), c[verts[j]][1]*n[1] + c[verts[j]][2]*n[2] + c[verts[j]][3]*n[3] - v[1]*n[1] - v[2]*n[2] - v[3]*n[3]]; 
end do; 
for i to nops(temp) do 
  if 0 < evalf(temp[i]) then 
    count1 := count1 + 1; 
  else 
    count2 := count2 + 1; 
  end if; 
end do; 
if count2 < count1 and count2 <> 0 or count1 < count2 and count1 <> 0 or count1 = count2 then 
  return false; 
else 
  return true; 
end if; 
end proc:


HasEdgeInPlane := proc(surface, face1, face2) 
local i, j, ed; ed := Edges(surface); 
for i to nops(face1) do 
  for j to nops(face2) do 
    if member([face1[i], face2[j]], ed) or member([face2[j], face1[i]], ed) then 
      return true; 
    end if; 
  end do; 
end do; 
return false; 
end proc:



TetraTorusTesting := proc(surface, face1, face2, coor, id)  #everthing
local m, n, o, p, q, r; 
o := HasVertexInPlane(surface, face1, face2, coor); 
p := HasEdgeInPlane(surface, face1, face2); 
if o = false and p = false then           #id=1: with mirror symmetry and without self intersection
  if id = 1 then 
    m := IsVertexfaithful(coor); 
    q := AreAllVerticesOnOneSideOfPlane(surface, face1, face2, coor); 
    n := HasSelfIntersections2(surface, coor, face1, face2,id); 
    if m = true and q = true and n = false then 
      return true; 
    end if; 
  elif id=2 then                              #id=2: with mirror symmetry with self intersection        
    m:=IsVertexfaithful(coor);        
    if m=true then          
      return true;        
    end if;     
  elif id=3 then                               #id=3: without mirror symmetry and without self intersection        
    r:=IsVertexfaithful_OnlyFaces(coor,face1,face2);        
    n:=HasSelfIntersections2(surface,coor,face1,face2,id);      
    if r=true and n=false then        
      return true;        
    end if;
  elif id = 4 then                             #id=4: without mirror symmetry with self intersection
    r := IsVertexfaithful_OnlyFaces(coor, face1, face2); 
    if r=true then 
      return true; 
    end if; 
  end if; 
end if; 
return false; 
end proc:





TetraTorusTest:=proc(surface,id)
local g,vertices,res,faces1,faces2,face1,face2,dets,solution,sol,aa,bb,i,j,coor,s,res1,bool;
vertices:=VerticesOfDegreeThree(surface);res:=[]; res1:=[];
faces1:=FacesOfVertex(surface,vertices[1]);
faces2:=FacesOfVertex(surface,vertices[2]);
for i in [1,2,3] do
  for j in [1,2,3] do
    face1:=faces1[i]; 
    face2:=faces2[j];
    if nops({op(face1)} intersect {op(face2)})=0 then 
      dets:=help_check(surface,face1,face2);
      solution := solve([dets[1] = 0, dets[2] = 0, dets[3] = 0],{a,b}); 
      for sol in solution do 
        aa:=subs(sol,a);bb:=subs(sol,b); 
        if  type(evalf(aa),float) and type(evalf(bb),float) then  
          coor:=evala(subs({a=aa,b=bb},CoordinateMatrix(surface,1,"listlist"=true)));  
          bool:=TetraTorusTesting(surface,face1,face2,coor,id);
          if aa<>0 and bb<>0 and bool=true then 
            #s:=NewSurface();DefineEmbedding(s,subs({a=aa,b=bb},CoordinateMatrix(surface,1,"listlist"=true)),"vertices"= Vertices(surface),"faces"=Faces(surface));
            print("res upated",face1,face2);
            res:=[op(res),sol]; res1:=[op(res1),[face1,face2]];       
          end if;
        end if; 
      end do; 
    end if;
  end do;
end do;
return res,res1; 
end proc:


findPossibleTetraTorus:=proc(cactiList,id)
local g,res,cac,s,i,sol,temp,solutions,aa,res2,bb,j,k,fac,facetemp,res3,temp1;res:=[];res2:=[];i:=1;  
for cac in cactiList do 
  facetemp := []; temp1 := []; temp:=[[cac]];  print("i=",i);j:=1;
  sol:=TetraTorusTest(cac,id); 
  for k from 1 to nops(sol[1]) do 
    aa:=subs(sol[1][k][1],a); print("j=",j);
    bb:=subs(sol[1][k][2],b); 
    if Im(aa)=0 and Im(bb)=0 then 
      temp1 := [op(temp1), sol[1][k]]; print("temp added");
      facetemp:=[op(facetemp),sol[2][k]]; print(facetemp);
    end if ;j:=j+1;
  end do;
  temp := [op(temp), temp1, facetemp]; 
  i:=i+1;
  if nops(temp[2])=0 then 
    res2:=[op(res2),cac]; 
  else 
    res:=[op(res),temp]; 
  end if ; 
end do ;
return res,res2;
end proc:


MirrorForTorus:=proc(surface,face1,face2)  
local i,j,v1,v2,v3,w1,w2,w3,vec1,vec2,normalvec,ww,kk,eqn_of_plane,verts,verts_of_faces1_2,vert_to_reflect,L,m,cc,newfac,maxvert,fac,f,surf,placeholder,reflected_verts,temp,jj,ij,n;  
cc:=CoordinateMatrix(surface,1,"listlist"=true);
v1 := cc[face1[1]];
v2 := cc[face1[2]];
v3 := cc[face1[3]];
w1 := cc[face2[1]];
w2 := cc[face2[2]];
w3 := cc[face12[3]];
vec1 := Vector(v2 - v1);
vec2 := Vector(w2 - v1);
normalvec := CrossProduct(vec1, vec2);
ww := x*normalvec[1] + y*normalvec[2] + z*normalvec[3] + k = 0;
kk := subs([x = v3[1], y = v3[2], z = v3[3]], ww);
solve(kk);
eqn_of_plane := subs(k = %, ww);
verts := Vertices(surface);
verts_of_faces1_2 := {op(face1), op(face2)};
vert_to_reflect := convert({op(verts)} minus verts_of_faces1_2, list);
for i to nops(vert_to_reflect) do
  L := <cc[vert_to_reflect[i]]> + l*normalvec;
  eval(eqn_of_plane, [x = L[1], y = L[2], z = L[3]]);
  m := solve(%, l);
  cc := [op(cc), convert(<cc[vert_to_reflect[i]]> + 2*m*normalvec, list)];
end do;
maxvert := max(verts);
n := nops(vert_to_reflect);
reflected_verts := [seq(j + maxvert, j = 1 .. n)];
fac := Faces(surface);
newfac := [];
for i from 1 to nops(fac) do   
  f:={op(fac[i])} minus {op(vert_to_reflect)};   
  if nops(f)<3 then      
    temp:={op(fac[i])} minus f;     
    placeholder:=[];     
    for j from 1 to nops(temp) do       
      for jj from 1 to nops(vert_to_reflect) do         
        if temp[j]=vert_to_reflect[jj] then           
          placeholder:=[op(placeholder),jj];         
        end if;        
      end do; 
    end do;      
    newfac:=[op(newfac),[op(f),op(map(ij->reflected_verts[ij],placeholder))]];    
  end if; 
end do;  
verts:=[op(verts),op(reflected_verts)];   
surf:=NewSurface();   
DefineEmbedding(surf,evala(cc),"faces"=[op(fac),op(newfac)],"vertices"=verts); 
return surf,evala(cc),Vertices(surf),Faces(surf); 
end proc:


torusfile_gen := proc(surfaces, id, filename) 
local i, j, cac, res, sol, temp, res2, fd, s, k, aa, bb, facetemp, l, cc, ss, temp1, refres, mir; 
i := 1; res := []; res2 := []; refres := []; 
fd := fopen(filename, WRITE); 
for s in surfaces do print("i" = i); 
  facetemp := []; temp1 := []; 
  cac := ConstructCactus(s); 
  temp := [[cac]]; j := 1; 
  sol := TetraTorusTest(cac, id); 
  for k to nops(sol[1]) do 
    aa := subs(sol[1][k][1], a); 
    bb := subs(sol[1][k][2], b); 
    if Im(aa) = 0 and Im(bb) = 0 then 
      temp1 := [op(temp1), sol[1][k]]; 
      facetemp := [op(facetemp), sol[2][k]]; 
    end if; print("j" = j); j := j + 1; 
  end do; 
  temp := [op(temp), temp1, facetemp]; 
  i := i + 1; 
  if nops(temp[2]) = 0 then 
    res2 := [op(res2), cac]; 
  else 
    fprintf(fd, "Number:"); fprintf(fd, String(i)); 
    for l to nops(temp[2]) do 
      ss := NewSurface(); 
      cc := evala(subs({a = subs(temp[2][l], a), b = subs(temp[2][l], b)}, CoordinateMatrix(op(temp[1]), 1, "listlist" = true))); 
      DefineEmbedding(ss, cc, "vertices" = Vertices(op(temp[1])), "faces" = Faces(op(temp[1]))); 
      mir := MirrorForTorus(ss, op(temp[3][l])); 
      fprintf(fd, "\n"); fprintf(fd, "Coordinates:"); fprintf(fd, String(mir[2])); 
      fprintf(fd, "\n"); fprintf(fd, "Vertices:"); fprintf(fd, String(mir[3])); 
      fprintf(fd, "\n"); fprintf(fd, "VerticesofAllFaces:"); fprintf(fd, String(mir[4])); 
      refres := [op(refres), mir[1]]; 
    end do; 
    res := [op(res), temp]; 
  end if; 
end do; 
fclose(fd); 
return refres; 
end proc:

TetraTorusTest_all_id:=proc(surface)
local g,vertices,res,faces1,faces2,face1,face2,dets,solution,sol,aa,bb,i,j,coor,s,res1,bool,k;
vertices:=VerticesOfDegreeThree(surface);res:=[]; res1:=[];
faces1:=FacesOfVertex(surface,vertices[1]);
faces2:=FacesOfVertex(surface,vertices[2]);
for i in [1,2,3] do
  for j in [1,2,3] do
    face1:=faces1[i]; 
    face2:=faces2[j];
    if nops({op(face1)} intersect {op(face2)})=0 then 
      dets:=help_check(surface,face1,face2);
      solution := solve([dets[1] = 0, dets[2] = 0, dets[3] = 0],{a,b}); 
      for sol in solution do 
        aa:=subs(sol,a);bb:=subs(sol,b); 
        if  type(evalf(aa),float) and type(evalf(bb),float) then  
          coor:=evala(subs({a=aa,b=bb},CoordinateMatrix(surface,1,"listlist"=true)));  
          for k in [1,2] do
            bool:=TetraTorusTesting(surface,face1,face2,coor,k);
            if aa<>0 and bb<>0 and bool=true then 
              #s:=NewSurface();DefineEmbedding(s,subs({a=aa,b=bb},CoordinateMatrix(surface,1,"listlist"=true)),"vertices"= Vertices(surface),"faces"=Faces(surface));
              print("res upated",face1,face2);
              res:=[op(res),sol]; res1:=[op(res1),[face1,face2]];       
            end if;
          end do;
        end if; 
      end do; 
    end if;
  end do;
end do;
return res,res1; 
end proc:


torusfile_gen_time := proc(surfaces, filename) 
local i, j, cac, res, sol, temp, res2, fd, s, k, aa, bb, facetemp, l, cc, ss, temp1, refres, mir; 
i := 1; res := []; res2 := []; refres := []; 
fd := fopen(filename, WRITE); 
for s in surfaces do print("i" = i); 
  facetemp := []; temp1 := []; 
  cac := ConstructCactus(s); 
  temp := [[cac]]; j := 1; 
  sol := TetraTorusTest_all_id(cac); 
  for k to nops(sol[1]) do 
    aa := subs(sol[1][k][1], a); 
    bb := subs(sol[1][k][2], b); 
    if Im(aa) = 0 and Im(bb) = 0 then 
      temp1 := [op(temp1), sol[1][k]]; 
      facetemp := [op(facetemp), sol[2][k]]; 
    end if; print("j" = j); j := j + 1; 
  end do; 
  temp := [op(temp), temp1, facetemp]; 
  i := i + 1; 
  if nops(temp[2]) = 0 then 
    res2 := [op(res2), cac]; 
  else 
    fprintf(fd, "Number:"); fprintf(fd, String(i)); 
    for l to nops(temp[2]) do 
      ss := NewSurface(); 
      cc := evala(subs({a = subs(temp[2][l], a), b = subs(temp[2][l], b)}, CoordinateMatrix(op(temp[1]), 1, "listlist" = true))); 
      DefineEmbedding(ss, cc, "vertices" = Vertices(op(temp[1])), "faces" = Faces(op(temp[1]))); 
      mir := MirrorForTorus(ss, op(temp[3][l])); 
      fprintf(fd, "\n"); fprintf(fd, "Coordinates:"); fprintf(fd, String(mir[2])); 
      fprintf(fd, "\n"); fprintf(fd, "Vertices:"); fprintf(fd, String(mir[3])); 
      fprintf(fd, "\n"); fprintf(fd, "VerticesofAllFaces:"); fprintf(fd, String(mir[4])); 
      refres := [op(refres), mir[1]]; 
    end do; 
    res := [op(res), temp]; 
  end if; 
end do; 
fclose(fd); 
return refres; 
end proc:


#torusfile_gen_time := proc(surfaces, filename) 

#read "/export3/home/tmp/maple_vani/Embeddings-of-wild-coloured-surfaces/Cacti10.g"
#sol := torusfile_gen(surfaces, 1, "test");
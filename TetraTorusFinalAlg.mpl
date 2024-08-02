restart:libname:=libname, "/home/data/archiv/daniel/maple/lib10": libname := libname, "/users/kirekod/maple/lib":
with(Involutive):with(LinearAlgebra):with(SimplicialSurfaceEmbeddings):



SimplicialSurfaceEmbeddingsOptions("radical"=false);

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
RemoveFace(surf,face1); RemoveFace(surf,face2);
#RemoveFace(surf,face1); RemoveFace(surf,face2);
return surf,evala(cc),Vertices(surf),Faces(surf); 
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
          coor:=evala(subs({a=aa,b=bb},CoordinateMatrix(surface,1,"listlist"=true,"radical"=false))); 
           
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


TetraTorusTesting := proc(surface, face1, face2, coor, id)  #everthing
local m, n, o, p, q, r; 
o := HasVertexInPlane(surface, face1, face2, coor); 
p := HasEdgeInPlane(surface, face1, face2); 
n:=HasSelfIntersections(surface,coor);
if o = false and p = false then           #id=1: with mirror symmetry and without self intersection
  if id = 1 then 
    m := IsVertexfaithful(coor); 
    n := HasSelfIntersections2(surface, coor, face1, face2,id); 
    if m = true and n = false then 
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




torusfile_new:=proc(surfaces,filename)
    local i,j,fd,cac,solution,k,co,verts,vof,res,s;
    res:=[];
    fd := fopen(filename, WRITE);
    cac:=map(j->ConstructCactus(j),surfaces);
    for id in [1,2,3,4] do
        for i from 1 to nops(cac) do
            solution:=findPossibleTetraTorus([cac[i]],id)[1];
            if nops(solution)<>0 then
                verts:= Vertices(op(solution[j][1])):
                vof:=Faces(op(solution[j][1]));
                resultCoordinates:=[];### 
                for k from 1 to nops(solution[j][2]) do
                    co:=subs({a = subs(solution[j][2][k], a), b = subs(solution[j][2][k], b)}, CoordinateMatrix(op(solution[j][1]), 1, "listlist" = true));
                    s:=NewSurface();
                    DefineEmbedding(s,co,"vertices"=verts,"faces"=vof);         
                    ####tests
                    ## all algebraic solutions are stored in tempRes 
                    ####
                end do;
                if nops(resultCoordinates )<>0 then 
                    fprintf(fd, "Simplicial Surface Number"); fprintf(fd, String(i));fprintf(fd, ":\n");
                    fprintf(fd, "vertices:="); fprintf(fd, String(verts)); fprintf(fd, ";\n");    
                    fprintf(fd, "\n"); fprintf(fd, "verticesoffaces:="); fprintf(fd, String(vof));fprintf(fd, ";\n");
                    for cc in resultCoordinates do  
                        fprintf(fd,"-------------------------------------\n");
                        fprintf(fd, "Coordinates:="); fprintf(fd, String(co));fprintf(fd, ";\n");
                    od; 
                    fprintf(fd,"-------------------------------------------------------------------------------------\n");
                end if;
            end if;
        end do;
    end do;
    fclose(fd);
    return res;
end proc:

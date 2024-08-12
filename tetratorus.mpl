#####################################################################################
###########
########### Beginning of Preamble
###########

#restart:
libname := libname, "/users/kirekod/maple/lib":
with(ListTools):with(Involutive):with(LinearAlgebra):with(SimplicialSurfaceEmbeddings):
with(combinat):
SimplicialSurfaceEmbeddingsOptions("radical",false);

#read "/Users/reymondakpanya/Desktop/github_repos/Embeddings-of-wild-coloured-surfaces/Cacti25.g":


###########
########### End of Preamble
###########
#####################################################################################

#####################################################################################
###########
########### Beginning of HelperFunctions
###########

help_norm := proc(v) 
  return sqrt(v[1]^2 + v[2]^2 + v[3]^2); 
end proc:

Wingtips:=proc(surface,e)
  local g,faces,f,butt,verts;
  faces:=Faces(surface);
  butt:=select(f->member(e[1],f) and member(e[2],f),faces);
  verts:=map(f->select(v->not(member(v,e)),f)[1],butt);
  return verts;
end proc:

Butterfly:=proc(surface,e)
  local g,wingtips,v;
  wingtips:=Wingtips(surface,e);
  return map(v->[e[1],e[2],v],wingtips);
end proc:

EqForFacesInSamePlane := proc(surface, face1, face2) 
  local c, v11, v12, v13, v21, v22, v23, d1, d2, d3;
  c := simplify(CoordinateMatrix(surface, 1, "listlist" = true, "radical"=false)); 
  v11 := c[face1[1]]; 
  v12 := c[face1[2]]; 
  v13 := c[face1[3]]; 
  v21 := c[face2[1]]; 
  v22 := c[face2[2]]; 
  v23 := c[face2[3]]; 
  d1 := Determinant(<<v12 - v11> | <v13 - v11> | <v21 - v11>>); 
  d2 := Determinant(<<v12 - v11> | <v13 - v11> | <v22 - v11>>); 
  d3 := Determinant(<<v12 - v11> | <v13 - v11> | <v23 - v11>>); 
  return [simplify(d1), simplify(d2), simplify(d3)]; 
end proc:

#AreFacesIdentified := proc(coor, face1, face2) # maybe delete and test
#  local g, i, c; 
#  for i in [1, 2, 3] do
#    if nops(selected(cc -> is(help_norm(cc - c[face1[i]]) =0), map(j -> c[face2[j]], [1, 2, 3]))) <> 1 then 
#      return false; 
#    end if; 
#  end do; 
#  return true; 
#end proc:

#position of x in l
pos:=proc(l,x)
  local g,i;
  for i from 1 to nops(l) do
    if l[i]=x then 
      return i;
    end if;
  od;
  return false;
end proc:
###########
########### End of Helperfunctions
###########
################################################################################################



################################################################################################
###########
########### Beginning of Testfunctions
###########

HasIdentifiedFaces:=proc(surface,face1,face2) #### test
  local g,coor,i,j,temp,f2; 
  coor:=CoordinateMatrix(surface,1,"listlist"=true); 
  f2:=[];
  for i from 1 to 3 do
    temp:=select(j->is(help_norm(coor[face1[i]]-coor[j])=0),face2); 
    if nops(temp)<>1 then 
      return false,[],[];
    end if:
    f2:=[op(f2),temp[1]];
  end do; 
  return true,face1,f2;
end proc:


HasSelfIntersections := proc(surface,coor) 
  local f, e,ee, v, v1, v2, w, w1, mat, sol, matin,ed,det;
  ed:=Edges(surface); 
  for f in Faces(surface) do 
    for e in select(ee -> not (ee[1] in f or ee[2] in f), ed) do 
      v := coor[f[1]]; 
      v1 := coor[f[2]] - coor[f[1]]; 
      v2 := coor[f[3]] - coor[f[1]]; 
      w := coor[e[1]]; 
      w1 := coor[e[2]] - coor[e[1]]; 
      mat := (<<v1> | <v2> | <-w1>>); 
      det:=verify(evala(Determinant(mat)),0,truefalse); print(evalf(Determinant(mat)),det);
      if det=false then 
        matin := MatrixInverse(mat); 
        sol:=(evala(matin . <w - v>)); 
          if is(sol[3]>=0) and is(sol[3]<=1) and is(sol[2]>=0) and is(sol[1]>=0) and is((sol[1] + sol[2]) <= 1) then 
            return true;
          end if;
      end if; 
    end do; 
  end do; 
  return false;
end proc:




HasSelfIntersections_old := proc(surface,coor) 
  local f, e, v, v1, v2, w, w1, mat, sol, matin,ed,sol_num,solu;
  #c:=CoordinateMatrix(surface,1,"listlist"=true,"radical"=false); #converting coor to radical might improve the algorithm
  ed:=Edges(surface); 
  for f in Faces(surface) do  
    for e in select(e -> not (e[1] in f or e[2] in f), ed) do 
      v := coor[f[1]]; 
      v1 := coor[f[2]] - coor[f[1]]; 
      v2 := coor[f[3]] - coor[f[1]]; 
      w := coor[e[1]]; 
      w1 := coor[e[2]] - coor[e[1]]; 
      mat := (<<v1> | <v2> | <-w1>>); #print(evalf(Determinant(mat)),verify(evala(Determinant(mat)),0,truefalse));
      if verify(evala(Determinant(mat)),0,truefalse)=false then #verify(evala(Determinant(mat)),0,truefalse)=false or evalf(Determinant(mat))<=10^(-7) then #is(evalf(Determinant(mat))<>0) then ###HERE!!
        matin := MatrixInverse(mat); 
        solu:=(matin . <w - v>);
        sol_num := evalf(solu); 
        #print("before sol", sol_num);
        if Im(sol_num[1])=0 and Im(sol_num[2])=0 and Im(sol_num[3])=0 then
          sol:=convert(solu,radical);
          if is(sol[3]>=0) and is(sol[3]<=1) and is(sol[2]>=0) and is(sol[1]>=0) and is((sol[1] + sol[2]) <= 1) then 
            #print("sol",sol_num);
            return true;
          end if;
        end if;
      end if; 
    end do; 
  end do; 
  return false;
end proc:


IsVertexfaithful := proc(coor)  ## test 
  local g, temp_num,temp_alg, i,c;
  #coor:=CoordinateMatrix(surface,1,"listlist"=true,"radical"=false);
  for c in coor do
    #temp_num:= select(v -> is(help_norm(evalf(c - v)) <= 10.0^(-8)), coor); 
     temp_alg:= select(v -> is(help_norm(evala(c - v))= 0), coor);
    #if 1 < nops(temp_num) then
      #temp_alg:= select(v -> is(help_norm(c - v)= 0), temp_num);
      if 1< nops(temp_alg) then 
        return false;
      end if;
    #end if; 
  end do;
  return true; 
end proc:


AreIsometricPolyhedra:=proc(s1, s2)
  local g,coor1,coor2,b1,b2,M1,M2,gram1,gram2,l,i;
  coor1:=CoordinateMatrix(s1,1,"listlist"=true, "radical"=false);l:=nops(coor1);
  coor2:=CoordinateMatrix(s2,1,"listlist"=true, "radical"=false);
  b1:=Barycenter(s1,1);
  b2:=Barycenter(s2,1);
  M1:=Matrix(map(i-><i-b1>,coor1));
  M2:=Matrix(map(i-><i-b1>,coor2));
  gram1:=Transpose(M1).M1;
  gram2:=Transpose(M2).M2;
  return Equal(gram1,gram2);
end proc:

IsometryRepresentatives:=proc(surfaceList) ## new 
  local g,res,s,ss,i,bool,l;
  res:=[];
  for s in surfaceList do
    bool:=true;
    i:=1;
    l:=nops(res);
    while bool and i< l+1 do
      ss:=res[i];
      if AreIsometricPolyhedra(s,ss) then
        bool:=false;
      end if;
      i:=i+1;
    end do;
    if bool then 
      res:=[op(res),s];
    end if;
  end do;
  return res;
end proc:

HasMirrorSymmetries:=proc(surface)
  local g,coor,edges,vec1,vec2,vec3,M,scalars,sc,bool,wings,faces,vertices,posScalars,negScalars,c,i,v,e;
  #print("starting mirror");
  coor:=simplify(CoordinateMatrix(surface,1,"vertices"=[$1..nops(Vertices(surface))],"listlist"=true, "radical"=false));
  edges:=Edges(surface);
  faces:=Faces(surface);
  vertices:=Vertices(surface);
  for e in edges do 
    wings:=Wingtips(surface,e);
    vec1:=<coor[e[2]]>-<coor[e[1]]>;
    vec2:=<(coor[wings[1]]+1/2*(coor[wings[2]]-coor[wings[1]]))-coor[e[1]]>;#vec1:=<coor[3]-coor[2]>;vec2:=<coor[5]-coor[2]>;
    vec3:=CrossProduct(vec1,vec2); 
    M:=Matrix([vec1,vec2,vec3]);  
    if is(Determinant(M)<>0) and is(VertexDegrees(surface)[e[1]] mod 2 =0) and is(VertexDegrees(surface)[e[2]] mod 2 =0) then
      M:=MatrixInverse(M); bool:=true;
      scalars:=map(c->M.<c-coor[e[1]]>,coor);posScalars:=select(c->is(c[3]>0),scalars);negScalars:=select(c->is(c[3]<0),scalars);
      if nops(negScalars) <> nops(posScalars) then 
        bool:=false; 
      end if;
      i:=1; 
      while bool and i<=nops(posScalars) do
        sc:=posScalars[i];
        if nops(select(c->is(c[1]=sc[1]) and is(c[2]=sc[2]) and is(c[3]=-sc[3]) ,negScalars))=0 then 
          bool:=false;
        else 
          i:=i+1; 
        end if; 
      end do;
      if bool then 
        return true; #e  
      end if;
    end if;
    vec1:=<coor[wings[1]]-(coor[e[1]]+1/2*(coor[e[2]]-coor[e[1]]))>; v:=<(coor[e[1]]+1/2*(coor[e[2]]-coor[e[1]]))>;
    vec2:=<coor[wings[2]]-(coor[e[1]]+1/2*(coor[e[2]]-coor[e[1]]))>;
    vec3:=CrossProduct(vec1,vec2);   
    M:=Matrix([vec1,vec2,vec3]);
    if Determinant(M)<>0 then
      M:=MatrixInverse(M); bool:=true;
      scalars:=map(c->M.(<c>-v),coor);posScalars:=select(c->is(c[3]>0),scalars);negScalars:=select(c->is(c[3]<0),scalars);
      if nops(negScalars) <> nops(posScalars) then 
        bool:=false; 
      end if;
      i:=1;
      while bool and i<=nops(posScalars) do
        sc:=posScalars[i];
        if nops(select(c->is(c[1]=sc[1]) and is(c[2]=sc[2]) and is(c[3]=-sc[3]) ,negScalars))=0 then 
          bool:=false;
        else i:=i+1;  end if;
      end do;
      if bool then 
        return true;#e  
      end if;
    end if;
  end do;
  return false;
end proc:

### tests whether an edge is incident to a number of edges that is not equal to 2 and 
## checks whether the euler characteristic is 0. ( neccessary but not sufficient !)

IsSimplicialTorus:=proc(surface)
  local g,edges,vertices,faces;
  edges:=Edges(surface);
  if nops(select(e->nops(FacesOfEdge(surface,e))=2,edges)) = nops(edges) then
    vertices:=Vertices(surface);
    faces:=Faces(surface);
    if nops(vertices)-nops(edges)+nops(faces)=0 then 
      return true;
    fi;
  fi; 
  return false;
end proc:

###########
########### End of Testfunctions
###########
################################################################################################

################################################################################################
###########
########### Beginning of combinatorics
###########

#DegreesOfVertices := proc(surface)
#  local g, i, f, degrees, count; 
#  degrees := []; 
#  for i in Vertices(surface) do 
#    count := 0; 
#    for f in Faces(surface) do 
#      if member(i, f) then 
#        count := count + 1; 
#      end if; 
#    end do; 
#    degrees := [op(degrees), count]; 
#  end do; 
#  return degrees; 
#end proc:

FacesOfVertex := proc(surface, v) 
  local g, res, faces, f; 
  res := []; 
  faces := Faces(surface); 
  for f in faces do 
    if f[1] = v or f[2] = v or f[3] = v then 
      res := [op(res), f]; 
    end if; 
  end do; 
  return res; 
end proc:

FacesOfEdge:=proc(surface,e)
  local g,fov1,fov2,f;
  fov1:=FacesOfVertex(surface,e[1]);
  fov2:=FacesOfVertex(surface,e[2]);
  return select(f->member(f,fov2),fov1);
end proc:
VerticesOfDegreeThree := proc(surface)
  local g, res, v,vertdeg,vert; 
  res := [];
  vert:=Vertices(surface);
  vertdeg:=VertexDegrees(surface);
  for v in vert do 
    if vertdeg[pos(vert,v)] = 3 then 
      res := [op(res), v]; 
    end if; 
  end do; 
  return res; 
end proc:


HasEdgeInPlane := proc(surface, face1, face2)
  local i, j, ed; 
  ed := Edges(surface); 
  for i to nops(face1) do 
    for j to nops(face2) do 
      if member([face1[i], face2[j]], ed) or member([face2[j], face1[i]], ed) then 
        return true; 
      end if; 
    end do; 
  end do; 
  return false; 
end proc:
###########
########### End of Combinatorics
###########
################################################################################################



################################################################################################
###########
########### Beginning of CactiConstruction
###########

AddTetra := proc(surf, i, vert, f)
  local g, c, v, v1, v2, v3, n, sol, vv,coor; 
  coor:=CoordinateMatrix(surf,i,"listlist"=true, "radical"=false); 
  v := <coor[f[1]]>;
  v1 := <coor[f[2]]> - <coor[f[1]]>;
  v2 := <coor[f[3]]> - <coor[f[1]]>;
  n := CrossProduct(v1, v2);
  v3 := <coor[vert]>; #print(v1,v2);
  vv := [t*n[1] + v3[1], t*n[2] + v3[2], t*n[3] + v3[3]];
  sol := solve(n[1]*vv[1] + n[2]*vv[2] + n[3]*vv[3] = v[1]*n[1] + v[2]*n[2] + v[3]*n[3], {t}); 
  return subs(sol, [2*t*n[1] + v3[1], 2*t*n[2] + v3[2], 2*t*n[3] + v3[3]]); 
  end proc:


AddT:=proc(surface,face,V)  
  local g,newV,vertices,c,vof,copyVOF,newVertex,s,l;
  RemoveFace(surface, convert(face, list));
  newV := max(Vertices(surface)) + 1;
  vertices := [op(Vertices(surface)), newV];
  c := CoordinateMatrix(surface, 1, "listlist" = true, "radical"=false);
  vof := Faces(surface);
  vof := {seq({vof[l][1], vof[l][2], vof[l][3]}, l = 1 .. nops(vof))};
  vof := (vof minus {face}) union {{newV, face[1], face[2]}, {newV, face[1], face[3]}, {newV, face[2], face[3]}};
  copyVOF := [seq([vof[l][1], vof[l][2], vof[l][3]], l = 1 .. nops(vof))];
  newVertex := AddTetra(surface, 1, V, [face[1], face[2], face[3]]);
  c := [op(c), newVertex];
  s := NewSurface();
  DefineEmbedding(s, c, "vertices" = vertices, "faces" = copyVOF);
  return s;
end proc:

# given our "symbol" cacsym for a cactus ( cacsym=[[face,vertex],...]) we want to construct a cacti 
# by iteratively reflecting the vertex 'vertex' through the plane spanned by the face 'face' 
ConstructCactus := proc(symbol) 
  local g, c, x, h, surf; 
  x := a; 
  h := b; 
  c := [[1/2, 0, 0], [-1/2, 0, 0], [(x^2 - 1)/(2 + 2*x^2), x/(1 + x^2), h], [-(x^2 - 1)/(2 + 2*x^2), -x/(1 + x^2), h]]; 
  surf := NewSurface();
  DefineEmbedding(surf, c, "faces" = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]], "vertices" = [1, 2, 3, 4]); 
  for g in symbol do 
    surf := AddT(surf, g[1], g[2]); 
  end do; 
  return surf; 
end proc:

###########
########### End of CactiConstruction
###########
################################################################################################

################################################################################################
###########
########### Beginning of TorusConstruction
###########
MirrorSimplicialSurface:=proc(surface,face1,face2)  ### formally Mirror for torus 
  local i,j,v1,v2,v3,w1,w2,w3,vec1,vec2,normalvec,ww,kk,eqn_of_plane,verts,verts_of_faces1_2,vert_to_reflect,L,m,cc,newfac,maxvert,fac,f,surf,placeholder,reflected_verts,temp,jj,ij,n;  
  cc:=CoordinateMatrix(surface,1,"listlist"=true, "radical"=false);
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
  DefineEmbedding(surf,cc,"faces"=[op(fac),op(newfac)],"vertices"=verts); 
  RemoveFace(surf,face1); RemoveFace(surf,face2);  
  if IsSimplicialTorus(surf) then 
    return surf;
  else 
    return surface;
  fi;
end proc:


IdentifyFaces:=proc(surface,face1,face2) ### identify faces ## new 
  local g,verts,i,j,newVerts,coor,faces,newFaces,f,temp,s,p; #print("checkcheck");
  verts:=Vertices(surface);
  newVerts:=[];
  for i from 1 to nops(verts) do ## identify vertices in face1 and face2
    for j from 1 to 3 do
      if verts[i]=face1[j] then
        verts[i]:=face2[j];
      end if;
    end do;
  end do; #print("newverts",newVerts,"verts",verts);
  for i from 1 to nops(verts) do ##set new vertices 
    if not(member(verts[i],newVerts)) then 
      newVerts:=[op(newVerts),verts[i]];
    end if;
  end do; #print("newverts",newVerts,"verts",verts);
  #newVerts:=[op({op(newVerts)})];
  coor:=CoordinateMatrix(surface,1,"vertices"=newVerts,"radical"=false);
  faces:=Faces(surface);
  newFaces:=[]; #print("f1",face1,"f2",face2);
  for f in faces do 
    if nops(select(i->member(i,face1),f))<>3 and nops(select(i->member(i,face2),f))<>3 then
      temp:=[];
      for i from 1 to 3 do
      p:=pos(face1,f[i]); #print("p",p);
        if p<>false then
          temp:=[op(temp),pos(newVerts,face2[p])];
#          temp:=[op(temp),face2[p]];
        else
          temp:=[op(temp),pos(newVerts,f[i])];
#          temp:=[op(temp),f[i]];
        end if;
      end do; #print("temp",temp,"f",f);
      newFaces:=[op(newFaces),[op({op(temp)})]];
    end if; 
  end do; #print("newFaces",newFaces);print("faces",faces);
  newFaces:=MakeUnique(newFaces); #print("updated newFaces",newFaces);
  newFaces:=select(f1->nops(select(f2->is(f1=f2),newFaces))=1,newFaces);#print("updated newFaces",newFaces);
  s:=NewSurface(); 
  DefineEmbedding(s,coor,"faces"=newFaces,"vertices"=[$1..nops(newVerts)]);  
  if IsSimplicialTorus(s) then 
    return s;
  else 
    return surface;
  fi;
end proc:

###########
########### End of TorusConstruction
###########
################################################################################################

ConstructTorusFromCactus:=proc(surface)
  local g,vertices,res,faces1,faces2,face1,face2,dets,solution,sol,aa,bb,i,j,coor,s,res1,bool,tempSurfaces,bool_mirror,bool_selfint,check_mirror,ss,coord,
  data_identify,f1,f2;
  vertices:=VerticesOfDegreeThree(surface); # exactly two vertices 
  faces1:=FacesOfVertex(surface,vertices[1]);
  faces2:=FacesOfVertex(surface,vertices[2]);
  tempSurfaces:=[[],[]];
  for i in [1,2,3] do
    for j in [1,2,3] do
      face1:=faces1[i]; 
      face2:=faces2[j]; 
      if nops({op(face1)} intersect {op(face2)})=0 and not(HasEdgeInPlane(surface, face1, face2)) then 
        dets:=EqForFacesInSamePlane(surface,face1,face2);
        solution := solve([dets[1] = 0, dets[2] = 0, dets[3] = 0],{a,b}); #print(nops([solution]),solution);
        #solution:=op(map(m->allvalues(solution[m],'implicit'),[$1..nops([solution])])); #print(nops([solution]),solution);
        for sol in solution do #print("check1");
          aa:=subs(sol,a);
          bb:=subs(sol,b); 
          if type(evalf(aa),float) and type(evalf(bb),float) and aa<>0 and bb<> 0 then 
            coor:=simplify(subs({a=aa,b=bb},CoordinateMatrix(surface,1,"listlist"=true,"radical"=false)));  
            s:=NewSurface();
            DefineEmbedding(s,coor,"faces"=Faces(surface),"vertices"= Vertices(surface));
            data_identify:=HasIdentifiedFaces(s,face1,face2); 
            if data_identify[1] then
              f1:=data_identify[2];
              f2:=data_identify[3]; 
              ss:=IdentifyFaces(s,f1,f2); 
              coord := simplify(CoordinateMatrix(ss,1,"listlist"=true, "radical"=false));
              check_mirror:=true;
            else
              ss:=MirrorSimplicialSurface(s,face1,face2); 
              coord := simplify(CoordinateMatrix(ss,1,"listlist"=true, "radical"=false));
              check_mirror:=false;
            end if;  
            if is(EulerCharacteristic(ss)=0) and IsVertexfaithful(coord) then 
                 #bool_selfint:=HasSelfIntersections(ss,coord); 
              if check_mirror then 
                bool_mirror:=HasMirrorSymmetries(ss); 
              else 
                bool_mirror:=true; # surface comes from imposing mirror symetries so nothing to be done 
              fi;
              if bool_mirror then # mirror
                tempSurfaces[1]:=[op(tempSurfaces[1]),ss]; 
              else  #non mirror
                tempSurfaces[2]:=[op(tempSurfaces[2]),ss];
              end if;  
            end if; 
          end if; 
        end do;  
      end if;
    end do;
  end do; 
  return tempSurfaces;
end proc:

ComputeCensus:=proc()
  local g,files,i;
  files:=["/export3/home/tmp/maple_vani/Embeddings-of-wild-coloured-surfaces/Cacti25.g","/export3/home/tmp/maple_vani/Embeddings-of-wild-coloured-surfaces/Cacti10.g"];
  for i from 1 to nops(files) do
    read files[i];
    torusfile(surfaces);
  end do:
  return true;
end proc:

#read "/export3/home/tmp/maple_vani/Embeddings-of-wild-coloured-surfaces/Cacti70.g":
#L := [[], [], [], [], [], []];
#
#
#for i from 1 to 80 do
#  cac:=ConstructCactus(surfaces[i]); print("i",i);
#  res:=ConstructTorusFromCactus(cac);
#  L[1]:=[op(L[1]),op(res[1])]; L[2]:=[op(L[2]),op(res[2])]; L[3]:=[op(L[3]),op(res[3])]; L[4]:=[op(L[4]),op(res[4])]; L[5]:=[op(L[5]),op(res[5])]; L[6]:=[op(L[6]),op(res[6])];
#end do:

################################################################################################
###########
########### Beginning of Functions to construct database
###########

torusfile:=proc(surfaces)
local i,j,cac,res,fd,tfile,l,L,nrFaces,co,filename;
filename:=["toriWithMirrorSymmetriesWith","toriWithoutMirrorSymmetriesWith"];
for i from 1 to nops(surfaces) do
    cac := ConstructCactus(surfaces[i]); print("i", i);
    L := ConstructTorusFromCactus(cac);    
    for j from 1 to 2 do
    	if nops(L[j])<>0 then
    		nrFaces:=nops(Faces(L[j][1]));
    		tfile:=cat(filename[j],convert(nrFaces,string),"Faces");
		fd := fopen(tfile, APPEND);
		fprintf(fd, "###################################################################################################################################### \n"); 
    		fprintf(fd, "SurfaceInfo:="); fprintf(fd, String(i)); fprintf(fd, ":\n");
    		fprintf(fd, "VerticesOfFaces:="); fprintf(fd, String(Faces(L[j][1]))); fprintf(fd, ";\n");
    		fprintf(fd, "Vertices:="); fprintf(fd, String(Vertices(L[j][1]))); fprintf(fd, ";\n");
    		for l from 1 to nops(L[j]) do
    			fprintf(fd, "##############################################################\n");
    			co:=evala(CoordinateMatrix(L[j][l],1,"listlist"=true,"radical"=true));
    			fprintf(fd, "coordinates_num:="); fprintf(fd, String(evalf(co))); fprintf(fd, ";\n");
    			fprintf(fd, "coordinates:="); fprintf(fd, String(co)); fprintf(fd, ";\n");
    		end do;
		fclose(fd);
	end if ;
    end do; 
end do;
end proc:


#torusfile_new:=proc(surfaces,filename)
#    local i,j,fd,cac,solution,k,co,verts,vof,res,s,id,resultCoordinates,cc;
#    res:=[];
#    fd := fopen(filename, WRITE);
#    cac:=map(j->ConstructCactus(j),surfaces);
#    for id in [1,2,3,4] do
#        for i from 1 to nops(cac) do
#            solution:=findPossibleTetraTorus([cac[i]],id)[1];
#            if nops(solution)<>0 then
#                verts:= Vertices(op(solution[j][1])):
#                vof:=Faces(op(solution[j][1]));
#                resultCoordinates:=[];### 
#                for k from 1 to nops(solution[j][2]) do
#                    co:=subs({a = subs(solution[j][2][k], a), b = subs(solution[j][2][k], b)}, CoordinateMatrix(op(solution[j][1]), 1, "listlist" = true, "radical"=false));
#                    s:=NewSurface();
#                    DefineEmbedding(s,co,"vertices"=verts,"faces"=vof);         
#                    ####tests
#                    ## all algebraic solutions are stored in tempRes 
#                    ####
#                end do;
#                if nops(resultCoordinates )<>0 then 
#                    fprintf(fd, "SurfaceInfo "); fprintf(fd, String(i));fprintf(fd, ":\n");
#                    fprintf(fd, "vertices:="); fprintf(fd, String(verts)); fprintf(fd, ";\n");    
#                    fprintf(fd, "\n"); fprintf(fd, "verticesoffaces:="); fprintf(fd, String(vof));fprintf(fd, ";\n");
#                    for cc in resultCoordinates do  
#                        fprintf(fd,"-------------------------------------\n");
#                        fprintf(fd, "Coordinates:="); fprintf(fd, String(co));fprintf(fd, ";\n");
#                    end do; 
#                    fprintf(fd,"-------------------------------------------------------------------------------------\n");
#                end if;
#            end if;
#        end do;
#    end do;
#    fclose(fd);
#    return res;
#end proc:
#
###########
########### End of Functions to construct database
###########
################################################################################################

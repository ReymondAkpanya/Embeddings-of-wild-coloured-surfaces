restart:libname:=libname, "/home/data/archiv/daniel/maple/lib10": libname := libname, "/users/kirekod/maple/lib":
with(Involutive):with(LinearAlgebra):with(SimplicialSurfaceEmbeddings):


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
  fclose(fd);
  return res;
end proc:


## script in gap
LoadPackage("simpl");;


myprint:= function(l)
    local i,t,m,temp;
    t:="[";
    temp:=Filtered([ 1 .. Length(l)],i->IsBound(l[i]));
    temp:=l{temp};
    for i in [1..Length(temp)] do
        t:=Concatenation(t,String(temp[i]));
        if i <> Length(temp) then 
            t:=Concatenation(t,",");
        else
       	    t:=Concatenation(t,"]");
        fi;	
    od;
    return t;
end;
f := Filename( DirectoryCurrent(), "deleteme.txt");
output := OutputTextFile( f, false );
if output = fail then
    Error(Concatenation("File ", String(file), " can't be opened.") );
fi;
SetPrintFormattingStatus( output, false );


s:=Octahedron();


#### TODO fancy stuff
str:=Concatenation("vertices:=",myprint(Vertices(s)),";\n");
str:=Concatenation(str,"edges:=",myprint(Edges(s)),";\n");
str:=Concatenation(str,"faces:=",myprint(Faces(s)),";\n");
str:=Concatenation(str,"vof:=",myprint(VerticesOfFaces(s)),";\n");


AppendTo(output,str);
CloseStream(output);
quit;


# goal: find all tetra helix with reflection symmetrie with up to n faces
#  1. construct all cacti with exactly two vertices of degree 3 with (n+1)/2 faces 
#  2. translate information from gap to maple (combinatorics -> embeddings)
#  3. Automize current maple code (solving equations for checking whether two faces lie in the same plane)
#  :) 

#given a list surfaceList with cacti (same number of faces, two vertices of degree 3) we want to construct the cacti with two more faces, two vertices of degree 3 (up to isomorphism)


###########################################################################
##########################################################################

FindCacti:=function(surfaceList,combList)
	local g,f,v,res,surf,verts,temp,neighbour,res2,tempSurfaces,orbFaces,faces,faceDeg,autGroup,vertices,help_IsomorphismRepresentatives,help_IsomorphismTest;
#################################################### helpfunction ################################################################################
#TODO analyze flame diagram of these function to see which function are expensive 
	help_IsomorphismTest:=function(surf,surfaceList)
		local g,vCounter,eCounter,fCounter,num_Geo,temp_surfaces,order_aut;
		## define invariants of surf
		vCounter:=ListCounter(CounterOfVertices(surf));
		eCounter:=ListCounter(CounterOfEdges(surf));
		fCounter:=ListCounter(CounterOfFaces(surf));
		#order_aut:=Order(AutomorphismGroup(surf));
		num_Geo:=Length(MaximalGeodesicPaths(surf));

		## use invariants to make the number of isomorphism tests smaller
		temp_surfaces:=Filtered(surfaceList,g->ListCounter(CounterOfVertices(g))=vCounter);
		temp_surfaces:=Filtered(temp_surfaces,g->ListCounter(CounterOfVertices(g))=vCounter);
		temp_surfaces:=Filtered(temp_surfaces,g->ListCounter(CounterOfVertices(g))=vCounter);
		#temp_surfaces:=Filtered(temp_surfaces,g->Order(AutomorphismGroup(g))=order_aut);
		temp_surfaces:=Filtered(temp_surfaces,g->num_Geo=Length(MaximalGeodesicPaths(g)));

		#num_wildcolourings:=Length(AllWildColouredSurfaces(surf)); #check how the other invariants perform and add them if neccessary 
		#num_tamecolourings:=Length(AllTameColouredSurfaces(surf));
		#temp_surfaces:=Filtered(temp_surfaces,g->num_wildcolourings=Length(AllWildColouredSurfaces(g)));
		#temp_surfaces:=Filtered(temp_surfaces,g->num_tamecolourings=Length(AllTameColouredSurfaces(g)));

		## isomorphism tests
		for g in temp_surfaces do 
			if IsIsomorphic(g,surf) then 
				return true;
			fi;
		od;
		return false;
	end;

	help_IsomorphismRepresentatives:=function(surfaceList,combList)
		local g,res,res2,n,i;
		res:=[];
		res2:=[];
		n:=Length(surfaceList);
		for i in [1..n] do
			if not help_IsomorphismTest(surfaceList[i],res) then 
				Add(res, surfaceList[i]);
				Add(res2,combList[i]);
			fi;
		od;
		return [res,res2];
	end;
#######################################################################################################################################################
	res:=[];
	res2:=[];
	for i in [1.. Length(surfaceList)] do
		surf:=surfaceList[i]; 
		vof:=VerticesOfFaces(surf);
		faceDeg:=FaceDegreesOfVertices(surf);
		vertices:=Filtered(Vertices(surf),v->faceDeg[v]=3);
		faces:=List(vertices,v->FacesOfVertex(surf,v));
		faces:=Union(faces); ## 6 elements
		autGroup:=AutomorphismGroupOnFaces(surf);
		orbFaces:=Orbits(autGroup,faces);
		faces:=List(orbFaces,g->g[1]);
		for f in faces do 
			Add(res,TetrahedralExtension(surf,f));
			verts:=vof[f];
			neighbour:=Intersection(List(verts,v->NeighbourVerticesOfVertexNC(surf,v)))[1];
			temp:=ShallowCopy(combList[i]);
			Add(temp,[vof[f],neighbour]);
			Add(res2,temp);
		od;
		
	od;
	#######
	return help_IsomorphismRepresentatives(res,res2);
end;

 ConstructCactus:=function(l)
s:=Tetrahedron();
 for g in l do
 s:=TetrahedralExtension(s,Position(VerticesOfFaces(s),g[1]));
 od;
 return s;
 end;


FindCacti_Loop:=function(n)
	local g,res,file;
	res:=[[T],[[]]];
	for i in [1..n] do 
		Print("i=",i,"\n");
		Print("surfaces=",Length(res[1]),"\n");
		res:=FindCacti(res[1],res[2]);
		file:="Cacti";
		file:=Concatenation(file,String(i));
		file:=Concatenation(file,".g");
		AppendTo(file,"surfaces:=");
		AppendTo(file,res[2]);
		AppendTo(file,";");
	od;
	return res;
end;




IsoRep:=function(surfList)
	local sList1,sList2,help_IsoList,help_IsomorphismTest,numGeo,fCounters,eCounters,vCounters;
	vCounters:=List(surfList,s->ListCounter(CounterOfVertices(s)));
	eCounters:=List(surfList,s->ListCounter(CounterOfEdges(s)));
	fCounters:=List(surfList,s->ListCounter(CounterOfFaces(s)));
	num_Geo:=List(surfList,s->Length(MaximalGeodesicPaths(s)));

	help_IsoList:=function(sList1,sList2)
		local g,res,n,i,surfaceList;
	
		res:=[];
		n:=Length(surfaceList);
	
		for i in [1..n] do
			if not help_IsomorphismTest(surfaceList[i],res) then 
				Add(res, surfaceList[i]);
			fi;
		od;
		return res;
	end;
end;


##

## new attempt 

writeCacti:=function()
	local g,i,iStream,str,edges,temp,res;
	res:=[];
	iStream:=InputTextFile("/home/data/akpanya/Konstruktion-und-Faltung-simplizialer-Flaechen/C24.g");
	for i in [1..58713] do
		str:=ReadLine(iStream);
		str:=SplitString(str,":")[2];
		str:=SplitString(str,"\n")[1];
		str:=SplitString(str,";");
		edges:=[];	

		for j in  [1..Length(str)] do
			temp:=SplitString(str[j],",");
			Add(edges,[Int(SplitString(temp[1]," ")[2]),Int(temp[2])]);
		od;


		## check degree
		degreeList:=Collected(List(Union(edges),i->Length(Filtered(edges,e->i in e))));
		if [3,2] in degreeList then 
			Add(res,edges);
		fi;
	od;
	Print(Length(res));
	vertices:=Union(res[1]);
	for i in [1..Length(res)] do
		all3Waists:=[];edges:=res[i];l:=Length(Union(edges));
		for i1 in [1.. l] do 
			for i2 in [i1+1..l] do
				for i3 in [i2+1..l] do 

					if IsSubset(edges,[[i1,i2],[i1,i3],[i2,i3]]) then 
						temp:=Filtered(vertices,i->Set([i,i1]) in edges and Set([i,i2]) in edges and Set([i,i3]) in edges );
						if Length(temp)=1 then 
						Add(all3Waists,[i1,i2,i3]);
					fi;fi;
				od;
			od;
		od;
		res[i]:=all3Waists;
	od;

	## construct vertices of faces


	return res;
end;

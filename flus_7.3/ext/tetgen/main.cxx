#include "tetgen.h" // Defined tetgenio, tetrahedralize().

extern "C" void tetgen(char tetgenopt[], int *nbnode, int *nbface, double bxyz[][3], int bf2n[][3], int bpatch[], int *nNode, int *nCell, int *nFace, double xyz[][3], int c2n[][4], int f2n[][3], int patch[] )
{
  tetgenio in, out, addin, bgmin;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int i;

  // All indices start from 1
  in.firstnumber = 1;

  // Boundary point informations 
  in.numberofpoints = *nbnode;
  in.pointlist = new REAL[in.numberofpoints * 3];
  for ( i = 0; i < in.numberofpoints; i++ ) {
	  in.pointlist[3*i]   = bxyz[i][0];
	  in.pointlist[3*i+1] = bxyz[i][1];
	  in.pointlist[3*i+2] = bxyz[i][2];
  }

  // Boundary facet informations
  in.numberoffacets = *nbface;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  for ( i = 0; i < in.numberoffacets; i++ ) {
	  f = &in.facetlist[i];
	  f->numberofpolygons = 1;
	  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	  f->numberofholes = 0;
	  f->holelist = NULL;
	  p = &f->polygonlist[0];
	  p->numberofvertices = 3;
	  p->vertexlist = new int[p->numberofvertices];
	  p->vertexlist[0] = bf2n[i][0];
	  p->vertexlist[1] = bf2n[i][1];
	  p->vertexlist[2] = bf2n[i][2];
	  in.facetmarkerlist[i] = bpatch[i];
  }

  // Tetrahedralize the PLC.
  tetrahedralize(tetgenopt, &in, &out, &addin, &bgmin);

  // Output New Grid Information
  *nNode = out.numberofpoints;
  *nFace = out.numberoftrifaces;
  *nCell = out.numberoftetrahedra;
  for ( i = 0; i < out.numberofpoints; i++ ) {
	  xyz[i][0] = out.pointlist[3*i];
	  xyz[i][1] = out.pointlist[3*i+1];
	  xyz[i][2] = out.pointlist[3*i+2];
  }
  for ( i = 0; i < out.numberoftetrahedra; i++ ) {
	  c2n[i][0] = out.tetrahedronlist[4*i];
	  c2n[i][1] = out.tetrahedronlist[4*i+1];
	  c2n[i][2] = out.tetrahedronlist[4*i+2];
	  c2n[i][3] = out.tetrahedronlist[4*i+3];
  }
  for ( i = 0; i < out.numberoftrifaces; i++ ) {
	  f2n[i][0] = out.trifacelist[3*i];
	  f2n[i][1] = out.trifacelist[3*i+1];
	  f2n[i][2] = out.trifacelist[3*i+2];
  }
  for ( i = 0; i < out.numberoftrifaces; i++ ) {
	  patch[i] = out.trifacemarkerlist[i];
  }

}

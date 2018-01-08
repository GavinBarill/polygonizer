#ifndef POLYGONIZER_HPP
#define POLYGONIZER_HPP
#include <iostream>

#include <string>
#include <cstdlib>
#include <vector>
/*
 * Adapted from 
 *
 * C code from the article
 * "An Implicit Surface Polygonizer"
 * http::www.unchainedgeometry.com/jbloom/papers/Polygonizer.pdf
 * by Jules Bloomenthal, jules@bloomenthal.com
 * in "Graphics Gems IV", Academic Press, 1994 */

/* implicit.c
 *     an implicit surface Polygonizer, translated from Mesa
 *     applications should call polygonize()
 *
 * To compile a test program for ASCII output:
 *     cc implicit.c -o implicit -lm
 *
 * To compile a test program for display on an SGI workstation:
 *     cc -DSGIGFX implicit.c -o implicit -lgl_s -lm
 *
 * Authored by Jules Bloomenthal, Xerox PARC.
 * Copyright (c) Xerox Corporation, 1991.  All rights reserved.
 * Permission is granted to reproduce, use and distribute this code for
 * any and all purposes, provided that this notice appears in all copies.  */

class Polygonizer
{
  public:
  // TYPES
  typedef struct point {		   /* a three-dimensional point */
      double x, y, z;		   /* its coordinates */
  } POINT;
  typedef struct test {		   /* test the function for a signed value */
      POINT p;			   /* location of test */
      double value;		   /* function value at p */
      int ok;			   /* if value is of correct sign */
  } TEST;
  typedef struct vertex {		   /* surface vertex */
      POINT position, normal;	   /* position and surface normal */
  } VERTEX;
  //typedef struct vertices {	   /* list of vertices in polygonization */
  //    int count, max;		   /* # vertices, max # allowed */
  //    VERTEX *ptr;		   /* dynamically allocated */
  //} VERTICES;
  typedef struct triangle {
          int i1, i2, i3;
  } TRIANGLE;
  //typedef struct triangles {
  //        int count, max;
  //        TRIANGLE *ptr;
  //} TRIANGLES;
  typedef struct corner {		   /* corner of a cube */
      int i, j, k;		   /* (i, j, k) is index within lattice */
      double x, y, z, value;	   /* location and function value */
  } CORNER;
  typedef struct cube {		   /* partitioning cell (cube) */
      int i, j, k;		   /* lattice location of cube */
      CORNER *corners[8];		   /* eight corners */
  } CUBE;
  typedef struct cubes {		   /* linked list of cubes acting as stack */
      CUBE cube;			   /* a single cube */
      struct cubes *next;		   /* remaining elements */
  } CUBES;
  typedef struct centerlist {	   /* list of cube locations */
      int i, j, k;		   /* cube location */
      struct centerlist *next;	   /* remaining elements */
  } CENTERLIST;
  typedef struct cornerlist {	   /* list of corners */
      int i, j, k;		   /* corner id */
      double value;		   /* corner value */
      struct cornerlist *next;	   /* remaining elements */
  } CORNERLIST;
  typedef struct edgelist {	   /* list of edges */
      int i1, j1, k1, i2, j2, k2;	   /* edge corner ids */
      int vid;			   /* vertex id */
      struct edgelist *next;	   /* remaining elements */
  } EDGELIST;
  typedef struct intlist {	   /* list of integers */
      int i;			   /* an integer */
      struct intlist *next;	   /* remaining elements */
  } INTLIST;
  typedef struct intlists {	   /* list of list of integers */
      INTLIST *list;		   /* a list of integers */
      struct intlists *next;	   /* remaining elements */
  } INTLISTS;
  typedef struct process {	   /* parameters, function, storage */
      /* implicit surface function */
      std::function< double(double,double,double) > function;
      /* triangle output function */
      std::function< bool(int,int,int,std::vector<VERTEX>&) > triproc;
      double size, delta;		   /* cube size, normal delta */
      int bounds;			   /* cube range within lattice */
      POINT start;		   /* start point on surface */
      CUBES *cubes;		   /* active cubes */
      //VERTICES vertices;		   /* surface vertices */
      std::vector<VERTEX> vertices;
      CENTERLIST **centers;	   /* cube center hash table */
      CORNERLIST **corners;	   /* corner value hash table */
      EDGELIST **edges;		   /* edge and vertex id hash table */
  } PROCESS;

  //
  // STATIC VARIABLES AND MACROS
  // These used to be macros... that might have been faster... I hope using
  // lambdas and ints does not change behavior...
  static const int HASHBIT = (5);
  static const int HASHSIZE = (size_t)(1<<(3*HASHBIT));   /* hash table size (32768) */
  static const int MASK = ((1<<HASHBIT)-1);
  static inline int HASH(const int i, const int j, const int k){return ((((((i)&MASK)<<HASHBIT)|((j)&MASK))<<HASHBIT)|((k)&MASK));};
  static inline int BIT(const int i, const int bit){return (((i)>>(bit))&1);};
  static inline int FLIP(const int i,const int bit){return ((i)^1<<(bit));}; /* flip the given bit of i */
  static const int RES = 10; /* # converge iterations    */
  static INTLISTS *cubetable[256];
  static const int TET = 0; /* use tetrahedral decomposition */
  static const int NOTET = 1; /* no tetrahedral decomposition  */
  static const int L = 0; /* left direction: = -x, -i  */
  static const int R = 1; /* right direction: = +x, +i  */
  static const int B = 2; /* bottom direction: -y, -j */
  static const int T = 3; /* top direction: = +y, +j  */
  static const int N = 4; /* near direction: = -z, -k  */
  static const int F = 5; /* far direction: = +z, +k  */
  static const int LBN = 0; /* left bottom near corner  */
  static const int LBF = 1; /* left bottom far corner   */
  static const int LTN = 2; /* left top near corner     */
  static const int LTF = 3; /* left top far corner      */
  static const int RBN = 4; /* right bottom near corner */
  static const int RBF = 5; /* right bottom far corner  */
  static const int RTN = 6; /* right top near corner    */
  static const int RTF = 7; /* right top far corner     */
  static const int LB = 0;  /* left bottom edge	*/
  static const int LT = 1;  /* left top edge	*/
  static const int LN = 2;  /* left near edge	*/
  static const int LF = 3;  /* left far edge	*/
  static const int RB = 4;  /* right bottom edge */
  static const int RT = 5;  /* right top edge	*/
  static const int RN = 6;  /* right near edge	*/
  static const int RF = 7;  /* right far edge	*/
  static const int BN = 8;  /* bottom near edge	*/
  static const int BF = 9;  /* bottom far edge	*/
  static const int TN = 10; /* top near edge	*/
  static const int TF = 11; /* top far edge	*/
  /*			edge: LB, LT, LN, LF, RB, RT, RN, RF, BN, BF, TN, TF */
  static constexpr int corner1[12]	   = {LBN,LTN,LBN,LBF,RBN,RTN,RBN,RBF,LBN,LBF,LTN,LTF};
  static constexpr int corner2[12]	   = {LBF,LTF,LTN,LTF,RBF,RTF,RTN,RTF,RBN,RBF,RTN,RTF};
  static constexpr int leftface[12]	   = {B,  L,  L,  F,  R,  T,  N,  R,  N,  B,  T,  F};
  /* face on left when going corner1 to corner2 */
  static constexpr int rightface[12]   = {L,  T,  N,  L,  B,  R,  R,  F,  B,  F,  N,  T};

  // DECLARATIONS
  static inline void converge(
    POINT * p1, 
    POINT * p2, 
    double v, 
    const std::function< double(double,double,double) > & function,
    POINT * p);
  static inline char *mycalloc(int nitems, int nbytes);

  /* docube: triangulate the cube directly, without decomposition */
  static int docube (CUBE * cube, PROCESS * p)
  {
      INTLISTS *polys;
      int i, index = 0;
      for (i = 0; i < 8; i++) if (cube->corners[i]->value > 0.0) index += (1<<i);
      for (polys = cubetable[index]; polys; polys = polys->next) {
          INTLIST *edges;
          int a = -1, b = -1, count = 0;
          for (edges = polys->list; edges; edges = edges->next) {
              CORNER *c1 = cube->corners[corner1[edges->i]];
              CORNER *c2 = cube->corners[corner2[edges->i]];
              int c = vertid(c1, c2, p);
              if (++count > 2 && ! p->triproc(a, b, c, p->vertices)) return 0;
              if (count < 3) a = b;
              b = c;
          }
      }
      return 1;
  }

  /* testface: given cube at lattice (i, j, k), and four corners of face,
  * if surface crosses face, compute other four corners of adjacent cube
  * and add new cube to cube stack */

  static inline void testface (
    int i,
    int j,
    int k,
    CUBE * old,
    int face,
    int c1,
    int c2,
    int c3,
    int c4,
    PROCESS * p);

  /* getedge: return vertex id for edge; return -1 if not set */
  
  static int getedge 
    (EDGELIST * table[], int i1, int j1, int k1, int i2, int j2, int k2)
  {
      EDGELIST *q;
      if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
  	int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
      };
      int index = HASH(i1, j1, k1)+HASH(i2, j2, k2) ;
      q = table[index];
      for (; q != NULL; q = q->next)
  	if (q->i1 == i1 && q->j1 == j1 && q->k1 == k1 &&
  	    q->i2 == i2 && q->j2 == j2 && q->k2 == k2)
  	    return q->vid;
      return -1;
  }
  /* setedge: set vertex id for edge */
  
  static inline void setedge (EDGELIST * table[], 
    int i1,
    int j1,
    int k1,
    int i2,
    int j2,
    int k2,
    int vid)
  {
      unsigned int index;
      EDGELIST *new_list;
      if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
  	int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
      }
      index = HASH(i1, j1, k1) + HASH(i2, j2, k2);
      new_list = (EDGELIST *) mycalloc(1, sizeof(EDGELIST));
      new_list->i1 = i1; new_list->j1 = j1; new_list->k1 = k1;
      new_list->i2 = i2; new_list->j2 = j2; new_list->k2 = k2;
      new_list->vid = vid;
      new_list->next = table[index];
      table[index] = new_list;
  }
  
  static inline void vnormal (POINT * point, PROCESS * p, POINT * v)
  {
      double f = p->function(point->x, point->y, point->z);
      v->x = p->function(point->x+p->delta, point->y, point->z)-f;
      v->y = p->function(point->x, point->y+p->delta, point->z)-f;
      v->z = p->function(point->x, point->y, point->z+p->delta)-f;
      f = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
      if (f != 0.0) {v->x /= f; v->y /= f; v->z /= f;}
  }

  /* vertid: return index for vertex on edge:
   * c1->value and c2->value are presumed of different sign
   * return saved index if any; else compute vertex and save */
  
  static int vertid (CORNER * c1,CORNER * c2, PROCESS * p)
  {
      VERTEX v;
      POINT a, b;
      int vid = getedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k);
      if (vid != -1) return vid;			     /* previously computed */
      a.x = c1->x; a.y = c1->y; a.z = c1->z;
      b.x = c2->x; b.y = c2->y; b.z = c2->z;
      converge(&a, &b, c1->value, p->function, &v.position); /* position */
      vnormal(&v.position, p, &v.normal);			   /* normal */
      p->vertices.push_back(v);
      vid = p->vertices.size()-1;
      setedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k, vid);
      return vid;
  }

  static int dotet (CUBE * cube, int c1, int c2, int c3, int c4, PROCESS * p)
  {
      CORNER *a = cube->corners[c1];
      CORNER *b = cube->corners[c2];
      CORNER *c = cube->corners[c3];
      CORNER *d = cube->corners[c4];
      int index = 0, apos, bpos, cpos, dpos, e1, e2, e3, e4, e5, e6;
      if ((apos = (a->value > 0.0))) index += 8;
      if ((bpos = (b->value > 0.0))) index += 4;
      if ((cpos = (c->value > 0.0))) index += 2;
      if ((dpos = (d->value > 0.0))) index += 1;
      /* index is now 4-bit number representing one of the 16 possible cases */
      if (apos != bpos) e1 = vertid(a, b, p);
      if (apos != cpos) e2 = vertid(a, c, p);
      if (apos != dpos) e3 = vertid(a, d, p);
      if (bpos != cpos) e4 = vertid(b, c, p);
      if (bpos != dpos) e5 = vertid(b, d, p);
      if (cpos != dpos) e6 = vertid(c, d, p);
      /* 14 productive tetrahedral cases (0000 and 1111 do not yield polygons */
          switch (index) {
          case 1:	 return p->triproc(e5, e6, e3, p->vertices);
          case 2:	 return p->triproc(e2, e6, e4, p->vertices);
          case 3:	 return p->triproc(e3, e5, e4, p->vertices) &&
                                          p->triproc(e3, e4, e2, p->vertices);
          case 4:	 return p->triproc(e1, e4, e5, p->vertices);
          case 5:	 return p->triproc(e3, e1, e4, p->vertices) &&
                                          p->triproc(e3, e4, e6, p->vertices);
          case 6:	 return p->triproc(e1, e2, e6, p->vertices) &&
                                          p->triproc(e1, e6, e5, p->vertices);
          case 7:	 return p->triproc(e1, e2, e3, p->vertices);
          case 8:	 return p->triproc(e1, e3, e2, p->vertices);
          case 9:	 return p->triproc(e1, e5, e6, p->vertices) &&
                                          p->triproc(e1, e6, e2, p->vertices);
          case 10: return p->triproc(e1, e3, e6, p->vertices) &&
                                          p->triproc(e1, e6, e4, p->vertices);
          case 11: return p->triproc(e1, e5, e4, p->vertices);
          case 12: return p->triproc(e3, e2, e4, p->vertices) &&
                                          p->triproc(e3, e4, e5, p->vertices);
          case 13: return p->triproc(e6, e2, e4, p->vertices);
          case 14: return p->triproc(e5, e3, e6, p->vertices);
          }
      return 1;
  }

  /* setcenter: set (i,j,k) entry of table[]
   * return 1 if already set; otherwise, set and return 0 */
  
  static int setcenter(CENTERLIST * table[], int i, int j, int k)
  {
      int index = HASH(i, j, k);
      CENTERLIST *newlist, *l, *q = table[index];
      for (l = q; l != NULL; l = l->next)
  	if (l->i == i && l->j == j && l->k == k) return 1;
      newlist = (CENTERLIST *) mycalloc(1, sizeof(CENTERLIST));
      newlist->i = i; newlist->j = j; newlist->k = k; newlist->next = q;
      table[index] = newlist;
      return 0;
  }

  /* setcorner: return corner with the given lattice location
    set (and cache) its function value */

  static CORNER *setcorner(PROCESS * p, int i, int j, int k)
  {
      /* for speed, do corner value caching here */
      CORNER *c = (CORNER *) mycalloc(1, sizeof(CORNER));
      int index = HASH(i, j, k);
      CORNERLIST *l = p->corners[index];
      c->i = i; c->x = p->start.x+((double)i-.5)*p->size;
      c->j = j; c->y = p->start.y+((double)j-.5)*p->size;
      c->k = k; c->z = p->start.z+((double)k-.5)*p->size;
      for (; l != NULL; l = l->next)
  	if (l->i == i && l->j == j && l->k == k) {
  	    c->value = l->value;
  	    return c;
  	    }
      l = (CORNERLIST *) mycalloc(1, sizeof(CORNERLIST));
      l->i = i; l->j = j; l->k = k;
      l->value = c->value = p->function(c->x, c->y, c->z);
      l->next = p->corners[index];
      p->corners[index] = l;
      return c;
  }



  static inline double RAND()
  {
    return ((rand()&32767)/32767.);
  }; /* random number between 0 and 1 */



  static inline TEST find (int sign, PROCESS * p, double x, double y, double z)
  {
      int i;
      TEST test;
      double range = p->size;
      test.ok = 1;
      for (i = 0; i < 10000; i++) {
          test.p.x = x+range*(RAND()-0.5);
          test.p.y = y+range*(RAND()-0.5);
          test.p.z = z+range*(RAND()-0.5);
          test.value = p->function(test.p.x, test.p.y, test.p.z);
          if (sign == (test.value > 0.0)) return test;
          range = range*1.0005; /* slowly expand search outwards */
      }
      test.ok = 0;
      return test;
  }

/* polygonize: polygonize the implicit surface function
 *   arguments are:
 *	 double function (x, y, z)
 *		 double x, y, z (an arbitrary 3D point)
 *	     the implicit surface function
 *	     return negative for inside, positive for outside
 *	 double size
 *	     width of the partitioning cube
 *	 int bounds
 *	     max. range of cubes (+/- on the three axes) from first cube
 *	 double x, y, z
 *	     coordinates of a starting point on or near the surface
 *	     may be defaulted to 0., 0., 0.
 *	 int triproc (i1, i2, i3, vertices)
 *		 int i1, i2, i3 (indices into the vertex array)
 *		 VERTICES vertices (the vertex array, indexed from 0)
 *	     called for each triangle
 *	     the triangle coordinates are (for i = i1, i2, i3):
 *		 vertices.ptr[i].position.x, .y, and .z
 *	     vertices are ccw when viewed from the out (positive) side
 *		 in a left-handed coordinate system
 *	     vertex normals point outwards
 *	     return 1 to continue, 0 to abort
 *	 int mode
 *	     TET: decompose cube and polygonize six tetrahedra
 *	     NOTET: polygonize cube directly
 *   returns error or NULL
 */
  // https://stackoverflow.com/a/6417182/148668
  template <typename LIST> 
  static void free_list(LIST * head)
  {
    LIST * tmp;
    while (head != NULL)
    {
       tmp = head;
       head = head->next;
       free(tmp);
    }
  };
  static inline std::string polygonize(
    const std::function< double(double,double,double) > & function,
    const double size,
    const int bounds,
    const double x,
    const double y,
    const double z,
    const std::function< bool(int,int,int,std::vector<VERTEX>&) > & triproc,
    int mode)
  {
    PROCESS p;
    int n, noabort;
    //std::function<CORNER*(PROCESS * p,int i, int j, int k)> setcorner;
    TEST in, out;
    //std::function<TEST(int sign,PROCESS * p,double x, double y, double z)> find;

    p.function = function;
    p.triproc = triproc;
    p.size = size;
    p.bounds = bounds;
    p.delta = size/(double)(RES*RES);
    /* allocate hash tables and build cube polygon table: */
    p.centers = new CENTERLIST*[HASHSIZE];
    p.corners = new CORNERLIST*[HASHSIZE];
    p.edges =	new EDGELIST*[2*HASHSIZE];
    //makecubetable();
    {
      int i, e, c, done[12], pos[8];
      for (i = 0; i < 256; i++) {
          for (e = 0; e < 12; e++) done[e] = 0;
          for (c = 0; c < 8; c++) pos[c] = BIT(i, c);
          for (e = 0; e < 12; e++)
              if (!done[e] && (pos[corner1[e]] != pos[corner2[e]])) {
                  INTLIST *ints = 0;
                  INTLISTS *lists = new INTLISTS[1];
                  int start = e, edge = e;
                  /* get face that is to right of edge from pos to neg corner: */
                  int face = pos[corner1[e]]? rightface[e] : leftface[e];
                  while (1) 
                  {
                    const auto & nextcwedge = [&](int edge, int face)->int
                    {
                      switch (edge) {
                        case LB: return (face == L)? LF : BN;
                        case LT: return (face == L)? LN : TF;
                        case LN: return (face == L)? LB : TN;
                        case LF: return (face == L)? LT : BF;
                        case RB: return (face == R)? RN : BF;
                        case RT: return (face == R)? RF : TN;
                        case RN: return (face == R)? RT : BN;
                        case RF: return (face == R)? RB : TF;
                        case BN: return (face == B)? RB : LN;
                        case BF: return (face == B)? LB : RF;
                        case TN: return (face == T)? LT : RN;
                        case TF: return (face == T)? RT : LF;
                      }
                      assert(false && "Should never get here");
                      return -1;
                    };

                      edge = nextcwedge(edge, face);
                      done[edge] = 1;
                      if (pos[corner1[edge]] != pos[corner2[edge]]) {
                          INTLIST *tmp = ints;
                          ints = new INTLIST[1];
                          ints->i = edge;
                          ints->next = tmp; /* add edge to head of list */
                          if (edge == start) break;
                          const auto & otherface = [&](const int edge, const int face)->int
                          {
                            int other = leftface[edge];
                            return face == other? rightface[edge] : other;
                          };
                          face = otherface(edge, face);
                      }
                  }
                  lists->list = ints; /* add ints to head of table entry */
                  lists->next = cubetable[i];
                  cubetable[i] = lists;
              }
      }
    }
    /* find point on surface, beginning search at (x, y, z): */
    srand(1);
    in = find(1, &p, x, y, z);
    out = find(0, &p, x, y, z);
    if (!in.ok || !out.ok) return "can't find starting point";
    converge(&in.p, &out.p, in.value, p.function, &p.start);

    /* push initial cube on stack: */
    p.cubes = new CUBES[1];
    p.cubes->cube.i = p.cubes->cube.j = p.cubes->cube.k = 0;
    p.cubes->next = NULL;

    /* set corners of initial cube: */
    for (n = 0; n < 8; n++)
	p.cubes->cube.corners[n] = setcorner(&p, BIT(n,2), BIT(n,1), BIT(n,0));

    p.vertices.clear();
    /* no vertices yet */

    setcenter(p.centers, 0, 0, 0);

    while (p.cubes != NULL) { /* process active cubes till none left */
        CUBE c;
        CUBES *temp = p.cubes;
        c = p.cubes->cube;

        noabort = mode == TET?
               /* either decompose into tetrahedra and polygonize: */
               dotet(&c, LBN, LTN, RBN, LBF, &p) &&
               dotet(&c, RTN, LTN, LBF, RBN, &p) &&
               dotet(&c, RTN, LTN, LTF, LBF, &p) &&
               dotet(&c, RTN, RBN, LBF, RBF, &p) &&
               dotet(&c, RTN, LBF, LTF, RBF, &p) &&
               dotet(&c, RTN, LTF, RTF, RBF, &p)
               :
               /* or polygonize the cube directly: */
               docube(&c, &p);
        if (! noabort) return "aborted";

        /* pop current cube from stack */
        p.cubes = p.cubes->next;
        // Of all the things to care about freeing from memory...
        free((char *) temp);
        /* test six face directions, maybe add to stack: */
        testface(c.i-1, c.j, c.k, &c, L, LBN, LBF, LTN, LTF, &p);
        testface(c.i+1, c.j, c.k, &c, R, RBN, RBF, RTN, RTF, &p);
        testface(c.i, c.j-1, c.k, &c, B, LBN, LBF, RBN, RBF, &p);
        testface(c.i, c.j+1, c.k, &c, T, LTN, LTF, RTN, RTF, &p);
        testface(c.i, c.j, c.k-1, &c, N, LBN, LTN, RBN, RTN, &p);
        testface(c.i, c.j, c.k+1, &c, F, LBF, LTF, RBF, RTF, &p);
    }

    // Shouldn't we free up all the memory we allocated?
    free_list(p.cubes);
    for(int i = 0;i<HASHSIZE;i++) free_list(p.centers[i]);
    delete[] p.centers;
    for(int i = 0;i<HASHSIZE;i++) free_list(p.corners[i]);
    delete[] p.corners;
    for(int i = 0;i<2*HASHSIZE;i++) free_list(p.edges[i]);
    delete[] p.edges;

    return "";
  }
};

// IMPLEMENTATIONS
/* mycalloc: return successful calloc or exit program */
inline char * Polygonizer::mycalloc(int nitems, int nbytes)
{
   char *ptr = (char*)calloc(nitems, nbytes);
   if (ptr != NULL) return ptr;
   assert(false);
   fprintf(stderr, "can't calloc %d bytes\n", nitems*nbytes);
   exit(1);
}
inline void Polygonizer::converge(
  POINT * p1, 
  POINT * p2, 
  double v, 
  const std::function< double(double,double,double) > & function,
  POINT * p)
{
  int i = 0;
  POINT pos, neg;
  if (v < 0) {
      pos.x = p2->x; pos.y = p2->y; pos.z = p2->z;
      neg.x = p1->x; neg.y = p1->y; neg.z = p1->z;
  }
  else {
      pos.x = p1->x; pos.y = p1->y; pos.z = p1->z;
      neg.x = p2->x; neg.y = p2->y; neg.z = p2->z;
  }
  while (1) {
      p->x = 0.5*(pos.x + neg.x);
      p->y = 0.5*(pos.y + neg.y);
      p->z = 0.5*(pos.z + neg.z);
      if (i++ == RES) return;
      if ((function(p->x, p->y, p->z)) > 0.0)
            {pos.x = p->x; pos.y = p->y; pos.z = p->z;}
      else {neg.x = p->x; neg.y = p->y; neg.z = p->z;}
  }
}

inline void Polygonizer::testface (
  int i,
  int j,
  int k,
  CUBE * old,
  int face,
  int c1,
  int c2,
  int c3,
  int c4,
  PROCESS * p)
{
    CUBE new_list;
    CUBES *oldcubes = p->cubes;
    static int facebit[6] = {2, 2, 1, 1, 0, 0};
    int n, pos = old->corners[c1]->value > 0.0 ? 1 : 0, bit = facebit[face];

    /* test if no surface crossing, cube out of bounds, or already visited: */
    if ((old->corners[c2]->value > 0) == pos &&
        (old->corners[c3]->value > 0) == pos &&
        (old->corners[c4]->value > 0) == pos) return;
    if (abs(i) > p->bounds || abs(j) > p->bounds || abs(k) > p->bounds) return;
    if (setcenter(p->centers, i, j, k)) return;

    /* create new_list cube: */
    new_list.i = i;
    new_list.j = j;
    new_list.k = k;
    for (n = 0; n < 8; n++) new_list.corners[n] = NULL;
    new_list.corners[FLIP(c1, bit)] = old->corners[c1];
    new_list.corners[FLIP(c2, bit)] = old->corners[c2];
    new_list.corners[FLIP(c3, bit)] = old->corners[c3];
    new_list.corners[FLIP(c4, bit)] = old->corners[c4];
    for (n = 0; n < 8; n++)
        if (new_list.corners[n] == NULL)
            new_list.corners[n] = setcorner(p, i+BIT(n,2), j+BIT(n,1), k+BIT(n,0));

    /*add cube to top of stack: */
    p->cubes = (CUBES *) mycalloc(1, sizeof(CUBES));
    p->cubes->cube = new_list;
    p->cubes->next = oldcubes;
}

// boo. global.
constexpr int Polygonizer::corner1[];
constexpr int Polygonizer::corner2[];
constexpr int Polygonizer::leftface[];
constexpr int Polygonizer::rightface[];
Polygonizer::INTLISTS *Polygonizer::cubetable[];


#endif

# Wonton Concepts      {#concepts}

## Wrappers 
Mesh and state or data associated with the mesh are fundamental
data structures in any physics based application. In general,
the application can have its own native data structures or 
use managers/frameworks for the mesh and associated state. 
Irrespectively, it becomes much harder to design generic
libraries targetting particular problems with the differences
in underlying mesh data structures and their access patterns. 
A key design aspect adopted by both Portage and Tangram is to
template on mesh and state wrapper types ensuring that any application
code can use them regardless of the underlying mesh data structures. 

The user is responsible for writing a mesh and state wrapper 
that provides basic information regarding how to access 
certain mesh and state components, e.g., how to get nodes, 
faces, cell data etc. Wonton is the library with such mesh
and state wrappers. Currently, Wonton contains wrappers to Jali mesh framework
and FleCSI Burton Mesh Specialization as well as a native
flat mesh/state wrappers (link to Dannys document). 

The mesh wrappers follow a Curiously Recurring Template Pattern 
CRTP pattern over a base class called AuxMeshTopology. The wrapper
implements methods for various mesh queries. It is derived from 
AuxMeshTopolgy that provides additional functionality to create 
side/wedge/corner types if they are not provided by the mesh framework. 
Since AuxMeshTopology needs to make some queries about the basic
topology, CRTP is used to allow this kind of mutual invocation from
the Derived wrapper type. 

The wrapper class must support cells, faces and nodes and adjacency
queries between these entities. In 2D, faces are the same as edges
and in 1D, faces are the same as nodes. The mesh class is expected
to support the following methods:

1. Basic 

    * int space\_dimension() const;  // dimensionality of mesh points (1, 2, 3)

    * int num\_owned\_cells() const;

    * int num\_ghost\_cells() const;

    * int num\_owned\_faces() const;

    * int num\_ghost\_faces() const;

    * int num\_owned\_nodes() const;

    * int num\_ghost\_nodes() const;

    * int get\_global\_id(int const id, Entity\_kind const kind) const;

2. Type information (see for details on Wonton types) 

    * Wonton::Entity\_type cell\_get\_type(int const cellid) const;

    * Wonton::Entity\_type node\_get\_type(int const nodeid) const;

    * Wonton::Element\_type cell\_get\_element\_type(int const cellid) const;

3. Topology queries
    * void cell\_get\_faces\_and\_dirs(int const cellid, std::vector<int> \*cfaces,
                              std::vector<int> *cfdirs) const;
      Here the cfdirs conveys the directions in which the faces are used by
      the cell. If the natural normal of the face points out of the cell, its
      direction should be returned as 1, if not, it should be returned as -1

    * void cell\_get\_nodes(int const cellid, std::vector<int> \*cnodes) const;

    * void face\_get\_nodes(int const faceid, std::vector<int> \*fnodes) const;

    * void face\_get\_cells(int const faceid, Wonton::Entity\_type etype,
                     std::vector<int> *fcells) const;

    * void node\_get\_cells(int const nodeid, Wonton::Entity\_type etype,
                     std::vector<int> *ncells) const;

    *  template<long D> void node\_get\_coordinates(int const nodeid, Wonton::Point<D> \*pp) const;

The AuxMeshTopology enhances the topological entities provided by the wrapper 
to build subcell entities - sides, corners and wedges. These non-standard
entities are required for remapping purposes when one or both meshes have elements
with non-planar faces. The definition of these entities are as follows:
* 1D: A side is a line segment from a node to the cell. Wedges and corners are the same
as the sides. 
* 2D: A side is a triangle formed by the two nodes of an edge/face and the cell center.
A wedge is half of a side formed by one node of the edge, the edge center and the 
cell center. A corner is a quadrilateral formed by the two wedges in a cell at a node. 
* 3D: A side is a tet formed by the two nodes of an edge, a face center and a cell center. 
A wedge is half a side, formed by a node of the edge, the edge center, the face 
center and the cell center. A corner is formed by all the wedges of a cell at a node.

At present, Wonton provides the following wrappers to external frameworks:
1. [Jali](http://github.com/lanl/jali)
2. [FleCSI Burton Specialization](http://github.com/laristra/flecsi-sp)

Wonton also provides two native mesh and state managers and their wrappers. 
1. Simple Mesh and State: A very light-weight, serial, 2D/3D Cartesian mesh and its state. 
3. Flat Mesh and State: (link to Dannys document)  


## Types
Wonton currently provides the following enumerated types:
1. _Entity\_kind_: The type of mesh entity. 
    * ALL\_KIND, ANY\_KIND, UNKNOWN\_KIND, NODE, EDGE, FACE, 
CELL, SIDE, WEDGE, CORNER, FACET, BOUNDARY\_FACE, PARTICLE  

2. _Entity\_type_: The parallel type of a given entity. 
    * TYPE\_UNKNOWN, DELETED, PARALLEL\_OWNED, PARALLEL\_GHOST, BOUNDARY\_GHOST, ALL

3. _Element\_type_: The cell type
    * UNKNOWN\_TOPOLOGY, TRI, QUAD, POLYGON, TET, PRISM, PYRAMID, HEX, POLYHEDRON

4. _Field\_type_: Whether the state i.e. the field is a mesh or multi-material field. 
    * UNKNOWN\_TYPE\_FIELD, MESH\_FIELD, MULTIMATERIAL\_FIELD

Wonton also provides types for Points, Vectors and Matrices. These are advanced types
with meaningfull operators for their particular type: 

* Point:
    - Operators: \[\], -p, p+t, c\*p, p/c
    - Methods: asV (convert Point to Vector) 
* Vector:
    - Operators: \[\], -V, V+U, V-U, sV, V/s, 
    - Methods: norm, one\_norm, max\_norm, normalize, is\_zero, dot, cross, MaxComponent 
* Matrix: 
    - Initialization: Using constant value, from vector of vectors, copy
    - Operators: M = M + C, M = M - C, M = c\*M, \[\], M = M\*_v_, M = A\*B
    - Methods: inverse, solve (QR decomposition based), data (return matrix data)  


## Algorithms

* SVD:  Double precision SVD decomposition routine. Takes an m-by-n matrix A and 
decomposes it into UDV, where U, V are left and right orthogonal transformation 
matrices, and D is a diagonal matrix of singular values. 

* lsfits: Computes a least squares gradient from a set of values. The first
  point is assumed to be the point where the gradient must be computed
  and the first value is assumed to the value at this reference point. 



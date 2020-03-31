# Simple Mesh Example    {#example}

Wonton provides a very crude mesh and state manager frameworks aptly
called `Simple_Mesh` and `Simple_State`.  The goal of these frameworks
is to show how one can wrap their favorite mesh and state manager for
use - _they should not be used in any production sense._

----

# Native Simple Mesh and State Manager

## Wonton::Simple_Mesh

This mesh framework is a non-distributed, 2D/3D, 
regular Cartesian mesh framework.  A Simple_Mesh is
constructed by specifying the extents of the box and the number of
cells in each direction.  The constructor then builds connectivity
information between cell, node, and face indices.  The ordering of
things like the nodes making up a cell, or the faces making up a cell
are specified consistently, but the choice of ordering does not
matter.  There are a few helper functions like
Simple_Mesh::cell_get_nodes() that will retrieve the indices
of some connected mesh entities given another mesh entity's index.

## Wonton::Simple_State

The state manager for Simple_Mesh is a collection
of named field data associated with a mesh entity type. For example, 
there can be a field named "density" which live on entities of the
type Wonton::CELL.  
The constructor for a Wonton::Simple_State takes a pointer to a
Wonton::Simple_Mesh so that it has access to things like the number
of nodes in the mesh.  Data are added to and retrieved from the state
manager via `add` and `get` methods.  Iterators are provided for the
map between `(name, location)` pairs and the data.

# Simple Mesh-State Wrappers

## Wonton::Simple_Mesh_Wrapper

This wraps the Wonton::Simple_Mesh framework.  Simple_Mesh, as its
name suggests, is quite simple and does not know about advanced mesh
entities or connectivities.  It lets [AuxMeshTopology](concepts.md###AuxMeshTopology) 
do the vast majority of the heavy lifting by automatically creating the advanced
mesh features.

Where possible, Simple_Mesh_Wrapper provides quick and efficient
answers to queries that AuxMeshTopology would otherwise compute 
more generally.  Two trivial examples are:

1. Wonton::Simple_Mesh_Wrapper::cell_get_type(), which determines the
   Wonton::Entity_type (e.g. PARALLEL_OWNED, PARALLEL_GHOST, etc.).
   We know Wonton::Simple_Mesh does not know anything about ghost
   information, so we simply return Wonton::Entity_type::PARALLEL_OWNED.
2. Wonton::Simple_Mesh_Wrapper::cell_get_element_type(), which
   determines the geometric shape of a given cell from one of the
   Wonton::Element_type's.  Simple mesh is only a 2d/3d, structured
   Cartesian mesh, so we return either Wonton::Element_type::QUAD 
   or  Wonton::Element_type::HEX.

There are a few other examples within Wonton::Simple_Mesh_Wrapper
where the AuxMeshTopology methods are overwritten to take advantage of
information that is cached within the Wonton::Simple_Mesh.  This is a
prime example of how the CRTP and AuxMeshTopology are intended to be
used.  

## Wonton::Simple_State_Wrapper

Wonton::Simple_State_Wrapper provides a set of 'add' and 'get' methods 
to retrieve data from the Simple_State as well as some error checking
for missing or duplicate data. 

# Other Meshes

Wonton provides other example meshes and wrappers to help users find or write
the tools they need.

## Adaptive Refinement Mesh and Wrapper

An extremely simplified mesh that mimics the structure of an adaptive mesh
refinement (AMR) grid.  This mesh should not be used in production, as it is
intended simply as a guide for the kinds of data structures and interfaces that
one might see in an AMR grid.  In addition, it only works in serial and not
across multiple processing elements.  The Adaptive Refinement Mesh is set up to
work with non-Cartesian coordinate systems, and can provide an example for
usage of various coordinate systems within Wonton and Portage.

An Adaptive Refinement Mesh requires three inputs:

1.  A refinement function (explained below).
2.  A point in space that defines the lowest corner of your mesh.
3.  A point in space that defines the highest corner of your mesh.

The outer bounds of the mesh are assumed to be an axis-aligned box extending
from the lowest corner to the highest corner.  Within that box, the mesh will
refine itself based on the refinement function.  It will evaluate the
refinement function at the center of each cell (initially the mesh is only a
single cell).  If the refinement function returns a value that is higher than
the current level of the cell (the initial cell is level zero), the cell is
refined by bisecting it along each axis.

The construction of a 2D Adaptive Refinement Mesh would look something like:

```
// coordinates extents
Point<2> lo = { 0.0, 0.0 };
Point<2> hi = { 1.0, 1.5 };
// the sizing function
auto refine = [&](Wonton::Point<2> const& p) {
  Point<2> norm;
  norm[0] = (p[0] - lo[0]) / (hi[0] - lo[0]);
  norm[1] = (p[1] - lo[1]) / (hi[1] - lo[1]);
  return static_cast<int>(std::ceil(2.0 + 3.0 * norm[0] + 1.0 * norm[1]));
};
// create the mesh
Wonton::Adaptive_Refinement_Mesh<2> mesh(refine, lo, hi);
```

This creates a 2D Adaptive Refinement Mesh that extends from (0.0,0.0) to
(1.0,1.5) and has a refinement function that creates finer cells at higher
values of x and y.

Technically this mesh only refines itself on setup, making it a
statically-refined mesh rather than an adaptively-refined mesh.  But it
demonstrates the sorts of interfaces and issues that an adaptively-refined mesh
needs to consider with respect to Portage.

## Direct Product Mesh and Wrapper

The Direct Product Mesh is a single-level, orthogonal mesh.  Along each axis
there is a set of points that specify the edges of the cells.  Taking a direct
product of these sets of points gives the set of corners for all cells in the
direct product mesh.  This allows the Direct Product Mesh to have variable cell
sizes along each axis, while still maintaining that simple, orthogonal
structure.  The Direct Product Mesh is set up to work with non-Cartesian
coordinate systems, and can provide an example for usage of various coordinate
systems within Wonton and Portage.

A Direct Product Mesh requires three inputs:

1.  A `std::array` of `std::vector` providing the list of edge points along
    each axis.
2.  An executor to define MPI communication details.
3.  The number of ghost cell layers.

The construction of a 2D Direct Product Mesh would look something like:

```
// instantiate an executor
Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);
// create the axis points
std::array<std::vector<double>,2> axis_points = {
    {0.0, 0.25, 0.75. 1.0},
    {0.0, 0.5, 1.5}
    };
// define how many ghost layers
int num_ghost_layers = 1;
// create the mesh
Wonton::Direct_Product_Mesh<1> mesh(axis_points, &executor, num_ghost_layers);
```

This creates a 2D Direct Product Mesh that is 3 cells by 2 cells, with MPI
communication, and a single layer of ghost cells for domain decomposition.

## Flat Mesh and Wrapper

## FleCSI Mesh and Wrapper

## Jali Mesh and Wrapper


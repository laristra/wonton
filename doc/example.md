# Simple Mesh Example    {#example}

Wonton provides a very crude mesh and state manager frameworks aptly
called `Simple_Mesh` and `Simple_State`.  The goal of these frameworks
is to show how one can wrap their favorite mesh and state manager for
use - _they should not be used in any production sense._

----

# The Mesh and State

## Wonton::Simple_Mesh

This mesh framework is a non-distributed (i.e. has no ghost
information), 2D/3D, regular Cartesian mesh framework.  A Simple_Mesh is
constructed by specifying the extents of the box and the number of
cells in each direction.  The constructor then builds connectivity
information between cell, node, and face indices.  The ordering of
things like the nodes making up a cell, or the faces making up a cell
are specified consistently, but the choice of ordering does not
matter.  There are a few helper functions like
Simple_Mesh::cell_get_nodes() that will retrieve the indices
of some connected mesh entities given another mesh entity's index.

## Wonton::Simple_State

The state manager for Simple_Mesh is essentially a collection
of field data specified by some name (e.g. "density") and location on
the Simple_Mesh where they live (e.g. Wonton::CELL).  The
constructor for a Wonton::Simple_State takes a pointer to a
Wonton::Simple_Mesh so that it has access to things like the number
of nodes in the mesh.  Data are added to and retrieved from the state
manager via `add` and `get` methods.  Iterators are provided for the
map between `(name, location)` pairs and the data.

## Wonton::Simple_Mesh_Wrapper

This wraps the Wonton::Simple_Mesh framework.  Simple_Mesh, as its
name suggests, is quite simple and does not know about advanced mesh
entities nor connectivities.  It lets AuxMeshTopology do the vast
majority of the heavy lifting by automatically creating the advanced
mesh features.

Where possible, Simple_Mesh_Wrapper provides quick and efficient
answers to queries that AuxMeshTopology would otherwise solve in a
general sense.  Two trivial examples are:

1. Wonton::Simple_Mesh_Wrapper::cell_get_type(), which determines the
   Wonton::Entity_type (e.g. PARALLEL_OWNED, PARALLEL_GHOST, etc.).
   We know Wonton::Simple_Mesh does not know anything about ghost
   information, so we simple return Wonton::Entity_type::PARALLEL_OWNED.
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

There is no equivalent of AuxMeshTopology for state wrappers.  This is
simply because the requirements of a state manager are much less
intense than those of the mesh framework.  In particular, the state
wrappers only need to know how to add data, get data, and query things
like the size of the data and where the data lives (Wonton::CELL or
Wonton::NODE).  Wonton::Simple_State_Wrapper exposes this
functionality and provides some error checking for missing or
duplicate data in the underlying Wonton::Simple_State object.


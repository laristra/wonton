/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_SUPPORT_WONTON_H_
#define WONTON_SUPPORT_WONTON_H_

// Autogenerated file that contains configuration specific defines
// like WONTON_ENABLE_MPI and WONTON_ENABLE_THRUST
#include "wonton-config.h"


#ifdef WONTON_ENABLE_THRUST

#include "thrust/device_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

#else  // no thrust

#include <boost/iterator/counting_iterator.hpp>

#include <vector>
#include <algorithm>
#include <string>

#endif

#ifdef WONTON_ENABLE_MPI
#include "mpi.h"
#endif

#ifdef WONTON_ENABLE_KOKKOS
#include <Kokkos_Macros.hpp>
#endif

/*
  @file wonton.h
  @brief Several utility types and functions within the Wonton namespace.
 */

/*
  @namespace Wonton
  The Wonton namespace houses all of the code within Wonton.

  Cells (aka zones/elements) are the highest dimension entities in a mesh
  Nodes (aka vertices) are lowest dimension entities in a mesh
  Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
  BOUNDARY_FACE is a special type of entity that is need so that process
  kernels can define composite vectors (see src/data_structures) on
  exterior boundary faces of the mesh only

  Wedges are special subcell entities that are a simplicial
  decomposition of cell. In 3D, a wedge is tetrahedron formed by one
  point of the edge, the midpoint of the edge, the "center" of the
  face and the "center" of the cell volume. In 2D, a wedge is a
  triangle formed by an end-point of the edge, the mid-point of the
  edge and the center of the cell. In 1D, wedges are lines, that are
  formed by the endpoint of the cell and the midpoint of the
  cell. There are two wedges associated with an edge of cell face in
  3D.

  Corners are also subcell entities that are associated uniquely with
  a node of a cell. Each corner is the union of all the wedges incident
  upon that node in the cell

  Facets are the boundary entity between two wedges in adjacent
  cells. In 3D, a facet is a triangular subface of the cell face
  shared by two wedges in adjacent cells. In 2D, a facet is half of
  an edge that is shared by two wedges in adjacent cells
 */
namespace Wonton {

typedef int64_t GID_t;  /*! Global ID type */



/// The type of mesh entity.
enum Entity_kind {
  ALL_KIND = -3,     /*!< All possible types */
  ANY_KIND = -2,     /*!< Any of the possible types */
  UNKNOWN_KIND = -1, /*!< Usually indicates an error */
  NODE = 0,
  EDGE,
  FACE,
  CELL,
  SIDE,
  WEDGE,
  CORNER,
  FACET,
  BOUNDARY_FACE,
  PARTICLE
};
constexpr int NUM_ENTITY_KIND = 13;

// Method to get Entity_kind in string form
inline std::string to_string(Entity_kind entkind) {
  static const std::string kind2string[NUM_ENTITY_KIND] =
    {"Entity_kind::ALL_KIND",
     "Entity_kind::ANY_KIND",
     "Entity_kind::UNKNOWN_KIND",
     "Entity_kind::NODE",
     "Entity_kind::EDGE",
     "Entity_kind::FACE",
     "Entity_kind::CELL",
     "Entity_kind::SIDE",
     "Entity_kind::WEDGE",
     "Entity_kind::CORNER",
     "Entity_kind::FACET",
     "Entity_kind::BOUNDARY_FACE",
     "Entity_kind::PARTICLE"};

  int itype = static_cast<int>(entkind)+3;
  return ((itype >= 0 && itype < NUM_ENTITY_KIND) ? kind2string[itype] :
          "INVALID Entity_kind");
}



// Parallel status of entity
/// The parallel type of a given entity.
enum Entity_type {
  TYPE_UNKNOWN = -1,
  DELETED = 0,
  PARALLEL_OWNED = 1,   /*!< Owned by this processor */
  PARALLEL_GHOST = 2,   /*!< Owned by another processor */
  BOUNDARY_GHOST = 3,   /*!< Ghost/Virtual entity on boundary */
  ALL  = 4              /*!< PARALLEL_OWNED + PARALLEL_GHOST + BOUNDARY_GHOST */
};
constexpr int NUM_ENTITY_TYPE = 6;

// Method to get Entity_type in string form
inline std::string to_string(Entity_type enttype) {
  static const std::string type2string[NUM_ENTITY_TYPE] =
    {"Entity_type::TYPE_UNKNOWN",
     "Entity_type::DELETED",
     "Entity_type::PARALLEL_OWNED",
     "Entity_type::PARALLEL_GHOST",
     "Entity_type::BOUNDARY_GHOST",
     "Entity_type::ALL"};

  int itype = static_cast<int>(enttype)+1;
  return ((itype >= 0 && itype < NUM_ENTITY_TYPE) ? type2string[itype] :
          "INVALID Entity_type");
}




/// Element (cell topology) type
enum Element_type {
  UNKNOWN_TOPOLOGY = 0,
  TRI,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHEDRON
};
constexpr int NUM_ELEMENT_TYPE = 9;

// Method to get Element_type in string form
inline std::string to_string(Element_type elemtype) {
  static const std::string type2string[NUM_ELEMENT_TYPE] =
    {"Element_type::UNKNOWN_TOPOLOGY",
     "Element_type::TRI",
     "Element_type::QUAD",
     "Element_type::POLYGON",
     "Element_type::TET",
     "Element_type::PRISM",
     "Element_type::PYRAMID",
     "Element_type::HEX",
     "Element_type::POLYHEDRON"};

  int itype = static_cast<int>(elemtype)+1;
  return ((itype >= 0 && itype < NUM_ELEMENT_TYPE) ? type2string[itype] :
          "INVALID Element_type");
}



/// Field type - whether it is mesh field or multi-material field
enum class Field_type {
  UNKNOWN_TYPE_FIELD = -1,
  MESH_FIELD,
  MULTIMATERIAL_FIELD
};
constexpr int NUM_FIELD_TYPE = 3;

inline std::string to_string(Field_type field_type) {
  static const std::string type2string[NUM_FIELD_TYPE] =
    {"Field_type::UNKNOWN_TYPE_FIELD",
     "Field_type::MESH_FIELD",
     "Field_type::MULTIMATERIAL_FIELD"};

  int itype = static_cast<int>(field_type)+1;
  return ((itype >= 0 && itype < NUM_FIELD_TYPE) ? type2string[itype] :
          "INVALID FIELD TYPE");
}


/// Layout of 2D input data to state - CELL_CENTRIC means the first
/// index is the cell and the second is the material; MATERIAL_CENTRIC
/// means the first index is the material and the second is the cell
enum class Data_layout {CELL_CENTRIC, MATERIAL_CENTRIC};
constexpr int NUM_DATA_LAYOUT = 2;

inline std::string to_string(Data_layout layout) {
  static const std::string type2string[NUM_DATA_LAYOUT] =
    {"Data_layout::CELL_CENTRIC", "Data_layout::MATERIAL_CENTRIC"};

  int itype = static_cast<int>(layout)+1;
  return ((itype >= 0 && itype < NUM_DATA_LAYOUT) ? type2string[itype] :
          "INVALID DATA LAYOUT");
}



/// Executor definition that lets us distinguish between serial and
/// parallel runs as well as furnish a communicator
struct Executor_type {
  virtual ~Executor_type() = default;
};

struct SerialExecutor_type : Executor_type {};  // for RTTI

#if !defined(WONTON_INLINE)
  #ifdef WONTON_ENABLE_KOKKOS
    #define WONTON_INLINE KOKKOS_INLINE_FUNCTION
  #else
    #define WONTON_INLINE inline
  #endif
#endif

#ifdef WONTON_ENABLE_MPI
struct MPIExecutor_type : Executor_type {
  explicit MPIExecutor_type(MPI_Comm comm) : mpicomm(comm) {}
  MPI_Comm mpicomm = MPI_COMM_WORLD;
};
#endif



#ifdef WONTON_ENABLE_THRUST

template<typename T>
    using vector = thrust::device_vector<T>;

template<typename T>
    using pointer = thrust::device_ptr<T>;

typedef thrust::counting_iterator<int> counting_iterator;
inline counting_iterator make_counting_iterator(int const i) {
  return thrust::make_counting_iterator(i);
}

template<typename InputIterator, typename OutputIterator,
         typename UnaryFunction>
inline OutputIterator transform(InputIterator first, InputIterator last,
                                OutputIterator result, UnaryFunction op) {
  struct thrust::execution_policy<thrust::system::omp::detail::tag> exec;
  return thrust::transform(exec, first, last, result, op);
}

template<typename InputIterator1, typename InputIterator2,
         typename OutputIterator, typename BinaryFunction>
inline OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                                InputIterator2 first2, OutputIterator result,
                                BinaryFunction op) {
  struct thrust::execution_policy<thrust::system::omp::detail::tag> exec;
  return thrust::transform(exec, first1, last1, first2, result, op);
}

template<typename InputIterator, typename UnaryFunction>
inline void for_each(InputIterator first, InputIterator last,
                              UnaryFunction f) {
  struct thrust::execution_policy<thrust::system::omp::detail::tag> exec;
  thrust::for_each(exec, first, last, f);
}

#else  // no thrust

template<typename T>
    using vector = std::vector<T>;

template<typename T>
    using pointer = T*;

typedef boost::counting_iterator<int> counting_iterator;
inline counting_iterator make_counting_iterator(int const i) {
  return boost::make_counting_iterator<int>(i);
}

template<typename InputIterator, typename OutputIterator,
    typename UnaryFunction>
inline OutputIterator transform(InputIterator first, InputIterator last,
                                OutputIterator result, UnaryFunction op) {
  return std::transform(first, last, result, op);
}

template<typename InputIterator1, typename InputIterator2,
         typename OutputIterator, typename BinaryFunction>
inline OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                                InputIterator2 first2, OutputIterator result,
                                BinaryFunction op) {
  return std::transform(first1, last1, first2, result, op);
}

template<typename InputIterator, typename UnaryFunction>
inline void for_each(InputIterator first, InputIterator last,
                     UnaryFunction f) {
  std::for_each(first, last, f);
}
#endif


struct Weights_t {
  Weights_t() : entityID(-1) {}
  Weights_t(int const entityID_in, std::vector<double> const& weights_in) :
      entityID(entityID_in), weights(weights_in) {}
  Weights_t(Weights_t const& source) :
      entityID(source.entityID), weights(source.weights) {}

  int entityID;
  std::vector<double> weights;
};

inline double pow2(double x) { return x*x; }


// Convert from some Portage types to MPI types
#ifdef WONTON_ENABLE_MPI
template<typename T> MPI_Datatype to_MPI_Datatype();
template<> inline MPI_Datatype to_MPI_Datatype<GID_t>() {
   return (sizeof(GID_t) == 8 ? MPI_LONG_LONG : MPI_INT);
}
#endif


}  // namespace Wonton

#endif  // WONTON_SUPPORT_WONTON_H_

/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef FLAT_MESH_WRAPPER_H_
#define FLAT_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include "wonton/mesh/AuxMeshTopology.h"
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

namespace Wonton {

using namespace Wonton;
/*!
  \class Flat_Mesh_Wrapper flat_mesh_wrapper.h
  \brief Flat_Mesh_Wrapper implements mesh methods

  Flat_Mesh_Wrapper stores mesh coordinates in a flat vector.
  It supports arbitrary polygons in 2D and arbitrary polyhedra
  in 3D.
*/

template <class T=double>
class Flat_Mesh_Wrapper : public AuxMeshTopology<Flat_Mesh_Wrapper<>> {
 public:

  //! Constructor
  Flat_Mesh_Wrapper() {};

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Flat_Mesh_Wrapper & operator=(Flat_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flat_Mesh_Wrapper() {};

  template<class Mesh_Wrapper>
  void initialize(Mesh_Wrapper& input)
  {

    dim_ = input.space_dimension();
    
    numOwnedCells_ = input.num_owned_cells();
    int numCells = numOwnedCells_ + input.num_ghost_cells();
    
    numOwnedNodes_ = input.num_owned_nodes();
    int numNodes = numOwnedNodes_ + input.num_ghost_nodes();   
     
    int numFaces = -1;
    if (dim_ == 3)
    {
      numOwnedFaces_ = input.num_owned_faces();
      numFaces = numOwnedFaces_ + input.num_ghost_faces();
    }
    
    // start clean
    reset_and_reserve(numCells, numFaces, numNodes);
    

    ///////////////////////////////////////////////////////
    // local copies we always need independent of dimension
    ///////////////////////////////////////////////////////
    
    // cell global ids, cell node counts, cell node listt
    for (unsigned int c=0; c<numCells; c++)
    {
      cellGlobalIds_.push_back(input.get_global_id(c, Entity_kind::CELL));

      std::vector<int> cellNodes;
      input.cell_get_nodes(c, &cellNodes);
      cellNodeCounts_.push_back(cellNodes.size());
      cellToNodeList_.insert(cellToNodeList_.end(),
                             cellNodes.begin(), cellNodes.end());
    }

    // node global ids
    for (unsigned int n=0; n<numNodes; ++n) {

      nodeGlobalIds_.push_back(input.get_global_id(n, Entity_kind::NODE));
    }
    

    ///////////////////////////////////////////////////////
    // local copies that are only required in 2D
    ///////////////////////////////////////////////////////
    
    if (dim_ ==2)
    {
    	// node coordinates
      for (unsigned int n=0; n<numNodes; ++n) {
        Wonton::Point<2> nodeCoord;
        input.node_get_coordinates(n, &nodeCoord);
        for (unsigned int j=0; j<2; ++j)
          nodeCoords_.push_back(nodeCoord[j]);
      }
    }
    
    
    ///////////////////////////////////////////////////////
    // local copies that are only required in 3D
    ///////////////////////////////////////////////////////
    
    if (dim_ == 3)
    {
    
      // cell face counts, cell face lists, cell face directions
      for (unsigned int c=0; c<numCells; ++c)
      {
        std::vector<int> cellFaces, cfDirs;
        input.cell_get_faces_and_dirs(c, &cellFaces, &cfDirs);
        cellFaceCounts_.push_back(cellFaces.size());
        cellToFaceList_.insert(cellToFaceList_.end(),
                               cellFaces.begin(), cellFaces.end());
        for (unsigned int j=0; j<cellFaces.size(); ++j)
          cellToFaceDirs_.push_back(cfDirs[j] >= 0);
      }

      // face global ids, face node counts, face node lists
      for (unsigned int f=0; f<numFaces; ++f)
      {
        faceGlobalIds_.push_back(input.get_global_id(f, Entity_kind::FACE));
        std::vector<int> faceNodes;
        input.face_get_nodes(f, &faceNodes);
        faceNodeCounts_.push_back(faceNodes.size());
        faceToNodeList_.insert(faceToNodeList_.end(),
                               faceNodes.begin(), faceNodes.end());
      }
      
      // node coordinates
      for (unsigned int n=0; n<numNodes; ++n) {
        Wonton::Point<3> nodeCoord;
        input.node_get_coordinates(n, &nodeCoord);
        for (unsigned int j=0; j<3; ++j)
          nodeCoords_.push_back(nodeCoord[j]);
      }
      
    } 
    
    // Complete the mesh. This step computes offsets and the inverse maps. 
    // We only need the compute_offsets routines before doing distribute, but 
    // it is possible (as in the tests) that someone will want to use a pre-
    // distribution flat mesh wrapper as a complete mesh, instead of the
    // underlying mesh, which would make more sense, so we will do this here
    // although it could be easily removed by adding 3 compute_offsets
    finish_init();

  }

  //! Finish mesh initialization, after initialize or MPI distribute
  void finish_init()
  {
    // Create all index maps
    make_index_maps();

    // Redo AuxMeshTopology info
    build_aux_entities();
  }

  // Create complete topologies for the various mesh entities.This topology is
  // defined by cells, faces and nodes. This is true in both dimensions 2 and 3
  // but obviously what is meant by these changes with dimension.
  // These topologies include the maps: cellToFace and faceToNode
  // and the skip map: cellToNode.
  // We also create the single inverse map going from nodeToCell.
  // This routine also computes offsets, which turn counts into partial sums so
  // that the relevant topologies, e.g. cell face indices, can be extracted from
  // the corresponding complete lists
  // Since faces are unique between cells, we also need to assign a direction to
  // each face in a cell
  void make_index_maps() {

    if (dim_ == 2)
    {
      // In 2D, we need to construct faces and their connectivity since the
      // boundary representation is just vertices around cells. A face in 2D is
      // and edge, and edges aren't specified, so we need to create them. Faces
      // also need to have directions and we need to make sure those are created
      // as well. 
      // In 2D we are given the skip topology cellToNode, and we need to construct
      // the adjacent topologies cellToFace and faceToNode.
      
      compute_offsets(cellNodeCounts_, &cellNodeOffsets_);

      // Clear and reserve
      cellToFaceList_.clear();
      cellToFaceList_.reserve(cellNodeCounts_.size());
      
      cellToFaceDirs_.clear();
      cellToFaceDirs_.reserve(cellNodeCounts_.size());
      
      faceToNodeList_.clear();
      faceToNodeList_.reserve(cellToNodeList_.size()); // slight underestimate
      
      // in 2D we directly specify the cell to node map. In order to complete
      // the topology we need to construct faces on the fly. Each face is 
      // defined by an (vertex, next_vertex) pair. There are as many faces(edges)
      // as vertices.
      
      // a temporary map that takes an edge node pair (v, next_v) and maps it to
      // a face id   
      std::map<std::pair<int, int>, int> nodePairToFaceTmp;
      
      // new face/edge counter
      int facecount = 0;
      
      for (unsigned int c=0; c<cellNodeCounts_.size(); ++c) {
        int offset = cellNodeOffsets_[c];
        int count = cellNodeCounts_[c];
        for (unsigned int i=0; i<count; ++i) {
          int n0 = cellToNodeList_[offset+i];
          int n1 = cellToNodeList_[offset+((i+1)%count)];
          // put nodes in canonical order
          int p0 = std::min(n0, n1);
          int p1 = std::max(n0, n1);
          auto npair = std::make_pair(p0, p1);
          auto it = nodePairToFaceTmp.find(npair);
          int face;
          if (it == nodePairToFaceTmp.end()) {
            // new face
            face = facecount;
            nodePairToFaceTmp.insert(std::make_pair(npair, face));
            faceToNodeList_.emplace_back(p0);
            faceToNodeList_.emplace_back(p1);
            ++facecount;
          }
          else {
            // existing face
            face = it->second;
          }
          cellToFaceList_.emplace_back(face);
          cellToFaceDirs_.push_back(n0 == p0);
        }  // for i

        if (c == numOwnedCells_ - 1) numOwnedFaces_ = facecount;

      }  // for c    
    }
    else if (dim_ == 3)
    {
      // In 3D, we are given the adjacent topologies cellToFace and faceToNode
      // and we only need to compute the skip topology cellToNode
      
      compute_offsets(faceNodeCounts_, &faceNodeOffsets_);
      compute_offsets(cellFaceCounts_, &cellFaceOffsets_);
      
      // Compute cell-to-node adjacency lists (3D only)
      cellNodeCounts_.clear();
      cellNodeCounts_.reserve(cellFaceCounts_.size());
      
      cellToNodeList_.clear();
      cellToNodeList_.reserve(cellFaceCounts_.size() * 4);
      
      for (unsigned int c=0; c<cellFaceCounts_.size(); ++c) {
        std::vector<int> cellfaces, dummy;
        cell_get_faces_and_dirs(c, &cellfaces, &dummy);
        std::set<int> cellnodes;

        for (auto const& f : cellfaces) {
          std::vector<int> facenodes;
          face_get_nodes(f, &facenodes);
          cellnodes.insert(facenodes.begin(), facenodes.end());
        }

        cellNodeCounts_.emplace_back(cellnodes.size());
        cellToNodeList_.insert(cellToNodeList_.end(),
                               cellnodes.begin(), cellnodes.end());
      }
      compute_offsets(cellNodeCounts_, &cellNodeOffsets_);
    }

    // Compute inverse node-to-cell adjacency lists
    int numNodes = nodeCoords_.size()/dim_;
    std::vector<std::set<int>> nodeToCellTmp(numNodes);
    for (unsigned int c=0; c<cellNodeCounts_.size(); ++c) {
      int offset = cellNodeOffsets_[c];
      int count = cellNodeCounts_[c];
      for (unsigned int i=0; i<count; ++i) {
        int n = cellToNodeList_[offset+i];
        nodeToCellTmp[n].insert(c);
      }
    }

    nodeCellCounts_.clear();
    nodeCellCounts_.reserve(numNodes);
    
    nodeToCellList_.clear();
    nodeToCellList_.reserve(cellToNodeList_.size());
    
    for (unsigned int n=0; n<numNodes; ++n) {
      const std::set<int>& nodes = nodeToCellTmp[n];
      nodeCellCounts_.emplace_back(nodes.size());
      nodeToCellList_.insert(nodeToCellList_.end(), nodes.begin(), nodes.end());
    }
    
    compute_offsets(nodeCellCounts_, &nodeCellOffsets_);

  } // make_index_maps

  //! Compute offsets from counts
  void compute_offsets(const std::vector<int>& counts,
                      std::vector<int>* offsets) {
    offsets->resize(counts.size());
    if (not counts.empty()) {
      (*offsets)[0] = 0;
      std::partial_sum(counts.begin(), counts.end()-1, offsets->begin()+1);
    }
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return numOwnedCells_;
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return (cellNodeCounts_.size() - numOwnedCells_);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return numOwnedNodes_;
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return (nodeCoords_.size()/dim_ - numOwnedNodes_);
  }

  //! Number of owned faces in the mesh
  int num_owned_faces() const {
    return numOwnedFaces_;
  }

  //! Number of ghost faces in the mesh
  int num_ghost_faces() const {
    if (dim_ == 2) {
      return (faceToNodeList_.size() / 2 - numOwnedFaces_);
    }
    else {
      return (faceNodeCounts_.size() - numOwnedFaces_);
    }
  }

  //! Coords of a node
  template <int D>
  void node_get_coordinates(int const nodeid, Point<D>* pp) const {
    assert(D == dim_);
    for (unsigned int j=0; j<dim_; j++)
      (*pp)[j] = nodeCoords_[nodeid*dim_+j];
  }

  //! Get the type of the cell - PARALLEL_OWNED or PARALLEL_GHOST
  Wonton::Entity_type cell_get_type(int const cellid) const {
    return (cellid < numOwnedCells_ ? PARALLEL_OWNED : PARALLEL_GHOST);
  }

  //! Get the type of the node - PARALLEL_OWNED or PARALLEL_GHOST
  Wonton::Entity_type node_get_type(int const nodeid) const {
    return (nodeid < numOwnedNodes_ ? PARALLEL_OWNED : PARALLEL_GHOST);
  }

  //! Get the element type of a cell - TRI, QUAD, POLYGON, TET, HEX,
  //! PRISM OR POLYHEDRON
  Wonton::Element_type cell_get_element_type(int const cellid) const {
    // FIXME
    return (dim_ == 2 ? POLYGON : POLYHEDRON);
  }

  //! Get cell faces and the directions in which they are used
  void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const {
    int offset, count;
    if (dim_ == 2)
    {
      offset = cellNodeOffsets_[cellid];
      count = cellNodeCounts_[cellid];
    }
    else
    {
      offset = cellFaceOffsets_[cellid];
      count = cellFaceCounts_[cellid];
    }
    cfaces->assign(&cellToFaceList_[offset],
                   &cellToFaceList_[offset+count]);
    cfdirs->clear();
    for (unsigned int j=0; j<count; ++j)
      cfdirs->push_back(cellToFaceDirs_[offset + j] ? 1 : -1);
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int const cellid, std::vector<int> *nodes) const {
    int offset = cellNodeOffsets_[cellid];
    int count = cellNodeCounts_[cellid];
    nodes->assign(&cellToNodeList_[offset],
                  &cellToNodeList_[offset+count]);
  }

  //! Get nodes of a face
  void face_get_nodes(int const faceid, std::vector<int> *fnodes) const {
    int offset, count;
    if (dim_ == 2)
    {
      offset = 2 * faceid;
      count = 2;
    }
    else
    {
      offset = faceNodeOffsets_[faceid];
      count = faceNodeCounts_[faceid];
    }
    fnodes->assign(&faceToNodeList_[offset],
                   &faceToNodeList_[offset+count]);
  }

  //! Get list of cells for a node
  void node_get_cells(int const nodeid, Entity_type const ptype,
                      std::vector<int> *cells) const {
    int offset = nodeCellOffsets_[nodeid];
    int count = nodeCellCounts_[nodeid];
    if (ptype == Entity_type::ALL)
      cells->assign(&nodeToCellList_[offset],
                    &nodeToCellList_[offset+count]);
    else {
      cells->clear();
      cells->reserve(count);
      for (int i = 0; i < count; i++) {
        int c = nodeToCellList_[offset+i];
        if (cell_get_type(c) == ptype)
          cells->push_back(c);
      }
    }
  }

  //! Coords of nodes of a cell
  template<int D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Wonton::Point<D>> *pplist) const {

    assert(D == dim_);
    std::vector<int> nodes;
    cell_get_nodes(cellid, &nodes);
    int cellNumNodes = nodes.size();
    pplist->resize(cellNumNodes);
    for (unsigned int i=0; i<cellNumNodes; ++i)
      node_get_coordinates(nodes[i], &((*pplist)[i]));
  }

  
  //! get coordinates
  std::vector<T>& get_coords() { return nodeCoords_; }

  /// cell --> node 
  //! get/set cell to node lists
  std::vector<int>& get_cell_to_node_list() { return cellToNodeList_; }
  
  void set_cell_to_node_list(std::vector<int>& cellToNodeList) 
  { cellToNodeList_ = cellToNodeList; }
  
  //! get/set cell node counts
  std::vector<int>& get_cell_node_counts() { return cellNodeCounts_; }
  
  void set_cell_node_counts(std::vector<int>& cellNodeCounts) 
  { cellNodeCounts_ = cellNodeCounts; }
  
  //! get cell node offsets
  std::vector<int>& get_cell_node_offsets() { return cellNodeOffsets_; }

  /// cell --> face
  //! get/set cell to face list
  std::vector<int>& get_cell_to_face_list() { return cellToFaceList_; }
  
  void set_cell_to_face_list(std::vector<int>& cellToFaceList) 
  { cellToFaceList_ = cellToFaceList; }

  //! get/set cell to face dirs
  std::vector<bool>& get_cell_to_face_dirs() { return cellToFaceDirs_; }
  
  void set_cell_to_face_dirs(std::vector<bool>& cellToFaceDirs) 
  { cellToFaceDirs_ = cellToFaceDirs; }

  //! get/set cell  face counts
  std::vector<int>& get_cell_face_counts() { return cellFaceCounts_; }

  void set_cell_face_counts(std::vector<int>& cellFaceCounts) 
  { cellFaceCounts_ = cellFaceCounts; }

  //! get cell face offsets
  std::vector<int>& get_cell_face_offsets() { return cellFaceOffsets_; }

  /// face --> node
  //! get/set face to node list
  std::vector<int>& get_face_to_node_list() { return faceToNodeList_; }
  
  void set_face_to_node_list(std::vector<int>& faceToNodeList) 
  { faceToNodeList_ = faceToNodeList; }

  //! get/set face node counts
  std::vector<int>& get_face_node_counts() { return faceNodeCounts_; }
 
  void set_face_node_counts(std::vector<int>& faceNodeCounts)
  { faceNodeCounts_ = faceNodeCounts; }

  //! get face node offsets
  std::vector<int>& get_face_node_offsets() { return faceNodeOffsets_; }

  /// node --> cell
  //! get node to cell list
  std::vector<int>& get_node_to_cell_list() { return nodeToCellList_; }

  //! get node cell counts
  std::vector<int>& get_node_cell_counts() { return nodeCellCounts_; }

  //! get node cell offsets
  std::vector<int>& get_node_cell_offsets() { return nodeCellOffsets_; }

 
  //! get global cell ids
  std::vector<GID_t>& get_global_cell_ids() { return cellGlobalIds_; }
  
  //! get global node ids
  std::vector<GID_t>& get_global_node_ids() { return nodeGlobalIds_; }
 
  //! get global face ids
  std::vector<GID_t>& get_global_face_ids() { return faceGlobalIds_; }
 
  void set_node_global_ids(std::vector<GID_t>& nodeGlobalIds) 
  { nodeGlobalIds_ = nodeGlobalIds; }

  //! set the number of owned cells
  void set_num_owned_cells(int numOwnedCells) 
  { numOwnedCells_ = numOwnedCells; }

  //! set the number of owned cells
  void set_num_owned_faces(int numOwnedFaces) 
  { numOwnedFaces_ = numOwnedFaces; }

  //! set the number of owned nodes
  void set_num_owned_nodes(int numOwnedNodes) 
  { numOwnedNodes_ = numOwnedNodes; }

  //! get spatial dimension
  int space_dimension() const { return dim_; }

  //! Get global ID of entities
  GID_t get_global_id(int ent, Entity_kind onwhat) const {
    switch (onwhat) {
      case Entity_kind::NODE: return nodeGlobalIds_[ent];
      case Entity_kind::FACE: return faceGlobalIds_[ent];
      case Entity_kind::CELL: return cellGlobalIds_[ent];
      default: return -1;
    }
  }

private:

  int dim_;
  int numOwnedCells_;
  int numOwnedFaces_;
  int numOwnedNodes_;

  std::vector<T>    nodeCoords_;
  std::vector<int>  cellToNodeList_;
  std::vector<int>  cellNodeCounts_;
  std::vector<int>  cellNodeOffsets_;
  std::vector<int>  cellToFaceList_;
  std::vector<bool> cellToFaceDirs_;  // unused in 2D (identical 
                                      // to cellNodeCounts_)
  std::vector<int>  cellFaceCounts_;  // unused in 2D (identical 
                                      // to cellNodeOffsets_)
  std::vector<int>  cellFaceOffsets_;
  std::vector<int>  faceToNodeList_;
  std::vector<int>  faceNodeCounts_;  // unused in 2D (always 2)
  std::vector<int>  faceNodeOffsets_; // unused in 2D (can be computed)
  std::vector<int>  nodeToCellList_;
  std::vector<int>  nodeCellCounts_;
  std::vector<int>  nodeCellOffsets_;
  
  std::vector<GID_t>  cellGlobalIds_;
  std::vector<GID_t>  faceGlobalIds_;
  std::vector<GID_t>  nodeGlobalIds_;
  

  void reset_and_reserve(int numCells, int numFaces, int numNodes){

    nodeCoords_.clear();
    nodeCoords_.reserve(numNodes*dim_);
    
    cellNodeCounts_.clear();
    cellNodeCounts_.reserve(numCells);
    
    cellToNodeList_.clear();
    cellToNodeList_.reserve(numCells*(dim_+1));
   
    cellGlobalIds_.clear();
    cellGlobalIds_.reserve(numCells);
    
    nodeGlobalIds_.clear();
    nodeGlobalIds_.reserve(numNodes);
    
    faceNodeCounts_.clear();
    faceToNodeList_.clear();
    faceGlobalIds_.clear();
    
    cellFaceCounts_.clear();
    cellToFaceList_.clear();
    cellToFaceDirs_.clear();
    
    // reserve (dim_+1) nodes per cell (lower bound)
    if (dim_ == 3)
    {
      // reserve know sizes
      faceGlobalIds_.reserve(numFaces);
      faceNodeCounts_.reserve(numFaces);
      cellFaceCounts_.reserve(numCells);
      
      // reserve lower bounds for sizes
      faceToNodeList_.reserve(numFaces*dim_);
      cellToFaceList_.reserve(numCells*(dim_+1));
      cellToFaceDirs_.reserve(numCells*(dim_+1));
    }
  }

}; // class Flat_Mesh_Wrapper

} // end namespace Wonton

#endif // FLAT_MESH_WRAPPER_H_

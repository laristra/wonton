/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef SIMPLE_STATE_MM_WRAPPER_H_
#define SIMPLE_STATE_MM_WRAPPER_H_

#include "wonton/state/state_manager.h"


/*!
  @file simple_state_wrapper.h
  @brief Definitions for a wrapper to Simple_State
*/
namespace Wonton {

using namespace Wonton;
  /*!
    @class Simple_State_Wrapper simple_state_wrapper.h
    @brief A thin wrapper that implements state methods for Simple_State needed
    by Wonton.
   */
template <class MeshWrapper>
class Simple_State_Wrapper: public StateManager<MeshWrapper> {

 public:
  /*!
    @brief Constructor for the state wrapper
   */
  Simple_State_Wrapper(const MeshWrapper& mesh, 
 			std::unordered_map<std::string,int> names={},
 			std::unordered_map<int,std::vector<int>> material_cells={}
 			) :StateManager<MeshWrapper>(mesh,names,material_cells) { }

  /// Assignment operator (disabled).
  Simple_State_Wrapper & operator=(Simple_State_Wrapper const &) = delete;

  /// Destructor.
  ~Simple_State_Wrapper() { }


 private:

};  // class Simple_State_Wrapper

}  // namespace Wonton 

#endif  // SIMPLE_STATE_MM_WRAPPER_H_

/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STATE_VECTOR_MULTI_H_
#define WONTON_STATE_VECTOR_MULTI_H_

#include <string>
#include <typeinfo>
#include <vector>
#include <unordered_map>

#include "wonton/support/wonton.h"
#include "wonton/state/state_vector_base.h"

namespace Wonton {  

/*!
	This class implements a state vector for a multi material field, meaning there
	are as many numbers per cell as materials in the cell. The current implementation
	is material dominant, meaning there is a list of cells for each material. 
	The shape of the data needs to be the same shape as the cell indices of 
	materials. Multi material fields can only be defined on cells.
*/
template <class T=double>
class StateVectorMulti : public StateVectorBase {

 public:
  
  /*!
    @brief Constructor for StateVectorMulti
    @param[in] data		The map of data. The key of the map is the material id.
    The value is the data defined on the material cells
  */
  StateVectorMulti(
      std::string name, 
      std::unordered_map<int, std::vector<T>> data = std::unordered_map<int, std::vector<T>>()
                   ) : StateVectorBase(name, Field_type::MULTIMATERIAL_FIELD, Entity_kind::CELL), 
                       data_(data) {}


  //! Destructor
  ~StateVectorMulti() {}

  // print
  std::ostream & print(std::ostream & os) const {
    os << "StateVectorMulti\n";
    return os;
  }

  // get the data type
  const std::type_info& data_type() const {
    const std::type_info& ti = typeid(T);
    return ti;
  }
	
	
  /*!
    @brief Return a reference to the data in the state vector.
    @return a reference to the vector of data in the state vector

    Return a reference to the data in the state vector.
  */
  std::unordered_map<int, std::vector<T>>& get_data() { return data_; }
		
		
  /*!
    @brief Return a const reference to the data in the state vector.
    @return a const reference to the vector of data in the state vector

    Return a const reference to the data in the state vector.
  */
  std::unordered_map<int, std::vector<T>> const & get_data() const { return data_; }
		
		
  /*!
    @brief Return a reference to the data in the state vector.
    @param[in] m	material id
    @return a reference to the vector of data in the state vector

    Return a reference to the data in the state vector. Note: this function
    has a potential side effect (which is a good thing, but users need to be
    aware). If the vector for key m does not exist, it will be created automatically.
    If you want to throw an error, use get_data().at(m) instead.
  */
  std::vector<T>& get_data(int m) { return data_[m]; }

  /*!
    @brief Return a const reference to the data in the state vector.
    @param[in] m	material id
    @return a const reference to the vector of data in the state vector

    Return a const reference to the data in the state vector. If you
    want to throw an error, use get_data().at(m) instead.
  */
  std::vector<T> const & get_data(int m) const { return data_.at(m); }

 private:
 
  std::unordered_map<int, std::vector<T>> data_;
 
};

}

#endif //WONTON_STATE_VECTOR_MULTI_H_


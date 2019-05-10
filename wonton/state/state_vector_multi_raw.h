/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STATE_VECTOR_MULTI_RAW_H_
#define WONTON_STATE_VECTOR_MULTI_RAW_H_

#include <string>
#include <typeinfo>
#include <memory>

#include "wonton/support/wonton.h"
#include "wonton/state/state_vector_base.h"

namespace Wonton {  

template <class T=double>
class StateVectorMultiRaw : public StateVectorBase {

 public:

  StateVectorMultiRaw(
      std::string name, 
      T** hdata=nullptr
                      ) : StateVectorBase(name, Field_type::MULTIMATERIAL_FIELD, Entity_kind::CELL), 
                          hdata_(hdata) {}


  //! Destructor
  ~StateVectorMultiRaw() {}

  // print
  std::ostream & print(std::ostream & os) const {
    os << "StateVectorMultiRaw\n";
    return os;
  }

  // get the data type
  const std::type_info& data_type() const {
    const std::type_info& ti = typeid(T);
    return ti;
  }
	
  /*!
    @brief Return a pointer to the data in the state vector.
    @return a pointer to the vector of data in the state vector

    Return a pointer to the data in the state vector.
  */
  T** get_data() { return hdata_; }

  /*!
    @brief Return a const pointer to the const data in the state vector.
    @return a pointer to the vector of data in the state vector

    Return a reference to the data in the state vector.
  */
  T const * const * get_data() const { return hdata_; }

 private:

  T** hdata_;

};

}

#endif //WONTON_STATE_VECTOR_MULTI_RAW_H_


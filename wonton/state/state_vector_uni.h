/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STATE_STATE_VECTOR_UNI_H_
#define WONTON_STATE_STATE_VECTOR_UNI_H_

#include <string>
#include <typeinfo>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/state/state_vector_base.h"

namespace Wonton {

/*!
        This class implements a state vector for a single material
        field, meaning there is only one number per mesh entity. The
        vector of data needs to be the same size as the number of
        entitites in the mesh. The field can be define on the cells,
        nodes, etc..
*/
template <class T = double>
class StateVectorUni : public StateVectorBase {
 public:
  /*!
    @brief Constructor for StateVectorUni
    @param[in] name         the name of the state vector
    @param[in] kind         the entity kind (CELL, NODE,...) of the data
    @param[in] data         the vector of data
  */
  StateVectorUni(
      std::string name,
      Entity_kind kind = Entity_kind::CELL,
      std::vector<T> data = std::vector<T>())
  : StateVectorBase(name, Field_type::MESH_FIELD, kind), data_(data) {}

  /*!
    @brief Constructor for StateVectorUni
    @param[in] name         the name of the state vector
    @param[in] kind         the entity kind (CELL, NODE,...) of the data
    @param[in] begin        pointer to the beginning of data to be copied
    @param[in] end          pointer to the end of the data to be copied

    Default arguments need to be at the end, so this signature is in a different
    order than the other constructor, and we need to be able to disambiguate
    with no default arguments.
  */
  StateVectorUni(
      std::string name,
      T const* begin,
      T const* end,
      Entity_kind kind = Entity_kind::CELL)
  : StateVectorBase(name, Field_type::MESH_FIELD, kind), data_(begin, end) {}


  //! Destructor
  ~StateVectorUni() {}

  // print
  std::ostream & print(std::ostream & os) const {
    os << "UniStateVector\n";
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
  std::vector<T>& get_data() { return data_; }

  /*!
    @brief Return a const reference to the data in the state vector.
    @return a const reference to the vector of data in the state vector

    Return a const reference to the data in the state vector.
  */
  std::vector<T> const & get_data() const { return data_; }

 private:
  std::vector<T> data_;
};  // class StateVectorUni

}  // namespace Wonton

#endif  // WONTON_STATE_STATE_VECTOR_UNI_H_

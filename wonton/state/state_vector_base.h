/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STATE_STATE_VECTOR_BASE_H_
#define WONTON_STATE_STATE_VECTOR_BASE_H_

#include <string>
#include <iostream>
#include <typeinfo>

#include "wonton/support/wonton.h"

namespace Wonton {

/*!
        This class implements a base class state vector. A
        StateVectorBase instance holds metadata common for all field
        types, such as the name of the field and on what mesh entity
        the field is defined. The primary reason for this base class
        is so that we can implement a single container that holds
        different conceptual field types such as single or multi
        material fields with different internal representations. We
        need a common base class pointer so we can store all the
        different types together. We can recast the pointer to the
        correct type whenever needed.
*/
class StateVectorBase {
 public:
  /*!
    @brief Constructor with a name
    @param name[in] Name of the StateVector
    @param type[in] Type of the StateVector (Field_type) (Single/Multi)
    @param kind[in] Kind of StateVector (CELL, NODE)
  */
  explicit StateVectorBase(std::string name, Field_type type,
                           Entity_kind kind = Entity_kind::CELL) :
  name_(name), type_(type), kind_(kind) {}


  //! Destructor

  virtual ~StateVectorBase() = default;

  //! Virtual methods

  virtual std::ostream & print(std::ostream & os) const {
    os << "Print not implemented for data type of StateVectorBase\n";
    return os;
  }


  virtual const std::type_info& data_type() const = 0;

  //! Query Metadata

  /*!
    @brief Return the name of the state vector.
    @return the name of the state vector

    Return the name of the state vector.
  */
  std::string get_name() const { return name_; }


  /*!
    @brief Return the field type [MESH_FIELD, MULTIMATERIAL_FIELD] of the state
    vector.
    @return the Field_type [MESH_FIELD, MULTIMATERIAL_FIELD] of the state vector

    Return the field type [MESH_FIELD, MULTIMATERIAL_FIELD] of the state vector.
  */
  Field_type get_type() const { return type_; }


  /*!
    @brief Return the entity kind [CELL, NODE] of the state vector.
    @return the Entity_kind [CELL, NODE] of the state vector

    Return the entity kind [CELL, NODE] of the state vector.
  */
  Entity_kind get_kind() const { return kind_; }


 protected:
  std::string name_;
  Field_type type_;
  Entity_kind kind_;
};

}  // namespace Wonton

#endif  // WONTON_STATE_STATE_VECTOR_BASE_H_

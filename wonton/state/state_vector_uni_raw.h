/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STATE_VECTOR_UNI_RAW_H_
#define WONTON_STATE_VECTOR_UNI_RAW_H_

#include <string>
#include <typeinfo>
#include <memory>

#include "wonton/support/wonton.h"
#include "wonton/state/state_vector_base.h"

namespace Wonton {  

template <class T=double>
class StateVectorUniRaw : public StateVectorBase {

	public:
  
		StateVectorUniRaw(
			std::string name, 
			Entity_kind kind=Entity_kind::CELL,
			T* pdata=nullptr
		) : StateVectorBase(name, Field_type::MESH_FIELD, kind),pdata_(pdata) {}


		//! Destructor
		~StateVectorUniRaw() {}

		// print
		std::ostream & print(std::ostream & os) const {
			os << "UniStateVector\n";
			return os;
		}
		
		// get the data type
		const std::type_info& data_type() {
			const std::type_info& ti = typeid(T);
			return ti;
		}
	
		/*!
			@brief Return a reference to the data in the state vector.
			@return a reference to the vector of data in the state vector

			Return a reference to the data in the state vector.
		*/
		T* get_data() { return pdata_; }

	private:
 
		T* pdata_;
 
};

}

#endif //WONTON_STATE_VECTOR_UNI_RAW_H_


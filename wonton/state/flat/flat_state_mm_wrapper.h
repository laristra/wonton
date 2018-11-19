/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/


#ifndef WONTON_STATE_FLAT_FLAT_STATE_MM_WRAPPER_H_
#define WONTON_STATE_FLAT_FLAT_STATE_MM_WRAPPER_H_

#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <type_traits>

#include "wonton/state/state_manager.h"
#include "wonton/support/wonton.h"
#include "wonton/state/state_vector_multi.h"

/*!
  @file flat_state_mm_wrapper.h
  @brief Definitions for a Flat State Wrapper that can do distributed MM remap
*/
namespace Wonton {

/*!
  @class Flat_State_Wrapper flat_state_wrapper.h
  @brief A state manager wrapper that allows redistribution of data across
  nodes.
*/
template <class MeshWrapper>
class Flat_State_Wrapper: public StateManager<MeshWrapper> {
 public:
  /*!
    @brief Constructor for the state wrapper
    @param[in] mesh           mesh wrapper
    @param[in] names          map from material names to material id
    @param[in] material_cells map from material id to vector of cells

    Constructor that takes the meshwrapper for the underlying mesh and
    two optional map arguments: the map from material name to id, and
    the map from material id to the cells containing that material
  */
  Flat_State_Wrapper(const MeshWrapper& mesh,
                     std::unordered_map<std::string, int> names = {},
                     std::unordered_map<int, std::vector<int>> material_cells =
                       {})
    : StateManager<MeshWrapper>(mesh, names, material_cells) { }


  /// Assignment operator (disabled).
  Flat_State_Wrapper & operator = (Flat_State_Wrapper const &) = delete;


  /// Destructor.
  ~Flat_State_Wrapper() { }


 /*!
  * @brief Initialize the state wrapper with another state wrapper and a list
  of names
  * @param[in] input another state wrapper, which need not be for Flat_State.
  * @param[in] var_names a list of state names to initialize
  *
  * Entities and sizes associated with the given name will be obtained from the
  input state wrapper.
  *
  * A name can be re-used with a different entity, but a name-entity
  combination must be unique.
  *
  * A name-entity combination must not introduce a new size for that entity if
  * it has previously been encountered.
  *
  * All existing internal data is forgotten.
  */
  template <class State_Wrapper>
  void initialize(State_Wrapper const & input,
                  const std::vector<std::string> var_names) {
                  
		// clear any existing state data                
		StateManager<MeshWrapper>::clear();  

	 	// get the number of materials
		int nmats = input.num_materials();
	
		// we may need to add volume fraction an centroid data to the list of 
		// distributed fields, so define a temp variable for these fields
		std::vector<std::string> distribute_var_names(var_names.begin(), var_names.end());
		
		// is this a multimaterial problem, if so, add material names and cell indices
		// to the flat state. Do this only once regardless of the number of fields.
		if (nmats>0){

			std::cout << "this is a multimaterial problem with " << nmats << " materials\n";
	
			// get the material names (for the following code to work, names() must
			// return the material names in a guaranteed order across processor nodes)
			std::vector<std::string> names = input.names();
	
			// create the name to id map			
			std::unordered_map<std::string, int> name_map;
			for (int i=0; i<nmats; ++i){
				name_map[names[i]]=i;
			}
	
			// add the names to the state manager
			StateManager<MeshWrapper>::add_material_names(name_map);

			// copy cell indices
			for (int m=0; m<nmats; ++m){		

				// get the cells for this material from the input state wrapper
				std::vector<int> mat_cells;
				input.mat_get_cells(m, &mat_cells);
				
				// add the cells to the flat state
				StateManager<MeshWrapper>::mat_add_cells(m, mat_cells);
				std::cout << "processing material " << m << " with " << mat_cells.size() << " cells (local ids)\n";
				for (const auto c : StateManager<MeshWrapper>::get_material_cells(m)) std:: cout << c << " ";
				std::cout << std::endl;  
			}  	
			
#ifdef HAVE_TANGRAM
			// if this is a multimaterial problem we will also need to distriube 
			// volume fraction and centroid
			distribute_var_names.emplace_back("mat_volfracs");
			distribute_var_names.emplace_back("mat_centroids");
#endif 

		}
	 
	 // add material fields to the flat state
   for (std::string varname : distribute_var_names) {
   
    // get the entity kind
    Entity_kind entity = input.get_entity(varname);   
   	
   // get the field type
    Field_type type = input.field_type(entity, varname);
    
    if (type==Wonton::Field_type::MULTIMATERIAL_FIELD){
    
    	std::cout << "****** got a MULTIMATERIAL_FIELD intialize: "<< varname <<std::endl;
    	
    	// We need a non-const state wrapper for the get_data_type function
    	State_Wrapper& input_nonconst = const_cast<State_Wrapper&>(input); 

    	// We need to create a state vector but the type is variable
    	const std::type_info& data_type = input_nonconst.get_data_type(varname);
    	
    	// construct the state vector type by if statement
    	// Is there a better/more idiomatic way to do this???
    	if (data_type == typeid(double)){
    	
    		std::cout << "found a double field\n";
    		
				// create a uni state vector
				std::shared_ptr<StateVectorMulti<double>> pv =
				    std::make_shared<StateVectorMulti<double>> (varname);
		      
				// add to database (inside loop, since pv is typed diffently below and 
				// pv has local scope
				StateManager<MeshWrapper>::add(pv);
				
		  	// loop over materials
		  	for (int m=0; m<nmats; ++m){
		  	     		
		  		// get the data from the input state wrapper
		  		double const* data;
		  		input.mat_get_celldata(varname, m, &data);
		  		
		  		// add the data to the flat state
		  		StateManager<MeshWrapper>::mat_add_celldata(varname, m, data);    		
		  		std::cout << varname << " data for material "<< m << std::endl;
		  		//auto pv = StateManager<MeshWrapper>::get(varname);
		  		//std::shared_ptr<StateVectorMulti<double>> pv1= std::static_pointer_cast<StateVectorMulti<double>>(pv);
		  		auto pv1 = StateManager<MeshWrapper>::get(varname);
		  		auto pv2 = std::dynamic_pointer_cast<StateVectorMulti<double>>(pv1);
		 		  for (const auto & d : pv2->get_data(m)) std:: cout << d << " ";
		   		std::cout << std::endl;
		  		
		  		
		  	}
    	} else if (data_type == typeid(Wonton::Point<2>)){
    	
    		std::cout << "found a Point<2>\n";
    		
				// create a uni state vector
				std::shared_ptr<StateVectorMulti<Wonton::Point<2>>> pv =
				    std::make_shared<StateVectorMulti<Wonton::Point<2>>> (varname);
		          		
				// add to database (inside loop, since pv is typed diffently below and 
				// pv has local scope
				StateManager<MeshWrapper>::add(pv);
				
		  	// loop over materials
		  	for (int m=0; m<nmats; ++m){
		  	     		
		  		// get the data from the input state wrapper
		  		Wonton::Point<2> const* data;
		  		input.mat_get_celldata(varname, m, &data);
		  		
		  		// add the data to the flat state
		  		StateManager<MeshWrapper>::mat_add_celldata(varname, m, data);    		
		  		std::cout << varname << " data for material "<< m << std::endl;
		  		//auto pv = StateManager<MeshWrapper>::get(varname);
		  		//std::shared_ptr<StateVectorMulti<double>> pv1= std::static_pointer_cast<StateVectorMulti<double>>(pv);
		  		auto pv1 = StateManager<MeshWrapper>::get(varname);
		  		auto pv2 = std::dynamic_pointer_cast<StateVectorMulti<Wonton::Point<2>>>(pv1);
		 		  for (const auto & d : pv2->get_data(m)) std:: cout << "(" << d[0] <<","<<d[1] << ") ";
		   		std::cout << std::endl;
		  		
		  		
		  	}
    	} else {
    		std::cout << "found a different field\n";
    	}
    	
    		
      std::cout << "****** end initialize\n";
    	
    } else {

		  // get pointer to data for state from input state wrapper
		  // DO WE NEED TO DO DATA INTROSPECTION TO FIND THE TYPE HERE (DOUBLE)
		  double const* data;
		  input.mesh_get_data(entity, varname, &data);

		  // Note, this get_data_size is of the input wrapper, which is defined
		  // not the base class state manager get_data_size which is not implemented
		  size_t dataSize = input.get_data_size(entity, varname);

		  // create a uni state vector
		  std::shared_ptr<StateVectorUni<double>> pv =
		      std::make_shared<StateVectorUni<double>> (varname, data,
		                                                data+dataSize, entity);
		  // add to database
		  StateManager<MeshWrapper>::add(pv);
    
		}
		
  }
  }


  /*!
   @brief Get the number of data vectors
   @return       The number of state vectors
  */
  size_t get_num_vectors() {
    return StateManager<MeshWrapper>::state_vectors_.size(); }


  /*!
   @brief Turn the field into a vector of doubles
   @param[in] field_name	The field name
   @return The serialized vector of doubles 
  
  */
  std::vector<double> serialize(std::string field_name) {
  
  	// get the state vector base class shared pointer
    std::shared_ptr<StateVectorBase> pv =
        StateManager<MeshWrapper>::get(field_name);
    
    // get the data type    
    const std::type_info& data_type = pv->data_type(); 
    
    // get the field type
    const Wonton:: Field_type field_type = pv->get_type();
    
    if (field_type == Wonton::Field_type::MESH_FIELD){

    	// this is a mesh field

  		if ( data_type == typeid(double)){

  			// field data is doubles
    		return std::static_pointer_cast<StateVectorUni<double>>(pv)->get_data();   		
    	}
    	
    } else if (field_type == Wonton::Field_type::MULTIMATERIAL_FIELD){

    	// this is a multimaterial field
    
  		if ( data_type == typeid(Wonton::Point<2>)){
  		
  			// field data is wonton 2d points
  			
  			// get the ordered keys
  			std::vector<int> const ids = get_material_ids();
  			
  			// get the raw mm data
    		std::unordered_map<int, std::vector<Wonton::Point<2>>> mm_data =  
    			std::static_pointer_cast<StateVectorMulti<Wonton::Point<2>>>(pv)->get_data();
    		
    		// define the result
    		std::vector<double> result;
    		
    		// loop over material ids in the mm state in sort order
    		for (int id : ids){
    			for (auto& d : mm_data.at(id)){
    				result.emplace_back(d[0]);
    				result.emplace_back(d[1]);
    			}
    		}
    		
    		return result;
    		
    	}  else	if (data_type == typeid(double)){
  		
  			// field data is doubles
  			
  			// get the ordered keys
  			std::vector<int> const ids = get_material_ids();
  			
  			// get the raw mm data
    		std::unordered_map<int, std::vector<double>> mm_data =  
    			std::static_pointer_cast<StateVectorMulti<double>>(pv)->get_data();
    		
    		// define the result
    		std::vector<double> result;
    		
    		// loop over material ids in the mm state in sort order
    		for (int id : ids){
    			for (auto& d : mm_data.at(id)){
    				result.emplace_back(d);
    			}
    		}
    		
    		return result;
    	}
    }
  }


   /*!
   @brief Turn a vector of doubles into the correct shape object
   @param[in] field_name	The field name
   @param[in] flat_data		The serialized striped vector of doubles to be unpacked
   @param[in] all_material_ids The striped material ids from all nodes
   @param[in] all_material_shapes The striped material shapes from all nodes
  
  */
  void deserialize(std::string field_name, const std::vector<double> & flat_data,
  	const std::vector<int> & all_material_ids, const std::vector<int> & all_material_shapes) {
  
  	// get the state vector base class shared pointer
    std::shared_ptr<StateVectorBase> pv =
        StateManager<MeshWrapper>::get(field_name);
    
    // get the data type    
    const std::type_info& data_type = pv->data_type(); 
    
    // get the field type
    const Wonton:: Field_type field_type = pv->get_type();
    
    if (field_type == Wonton::Field_type::MESH_FIELD){
    
    	// this is a mesh field
    
  		if ( data_type == typeid(double)){
  		
  			// copy the mesh data into the state vector
    		std::static_pointer_cast<StateVectorUni<double>>(pv)->get_data()=flat_data;
    	}
    	
    } else if (field_type == Wonton::Field_type::MULTIMATERIAL_FIELD){
    
    	// this is a multimaterial field
    
  		if ( data_type == typeid(Wonton::Point<2>)){
  		
  			// field data is 2d wonton points
  			
  			// allocate the vector data
  			std::vector<Wonton::Point<2>> correctly_typed_data(flat_data.size()/2);
  			
  			// convert to the correct type
  			for (int i=0; i<flat_data.size(); i+=2){
  				correctly_typed_data[i/2] = Wonton::Point<2>(flat_data[i], flat_data[i+1]);
  			}
  			
  			// allocate the data for unpacking
	      std::unordered_map<int,std::vector<Wonton::Point<2>>> material_data;
	      
    	
				/////////////////////////////////////////////////////////
		    // We need to turn the flattened material cells into a correctly shaped
		    // ragged right structure for use as the material cells in the flat
		    // state wrapper. Just as in the flat state mesh field (and associated
		    // cell ids), we aren't removing duplicates, just concatnating by material
		    /////////////////////////////////////////////////////////
		    
		    
		    // reset the running counter
		    int running_counter=0;
		    
		    // loop over material ids on different nodes
		    for (int i=0; i<all_material_ids.size(); ++i){
		    
		    	// get the current working material
		    	int mat_id = all_material_ids[i];
		    	
		    	// get the current number of material cells for this material
		    	int nmat_cells = all_material_shapes[i];
		    	
		    	// get or create a reference to the correct material cell vector
		    	auto& these_material_data = material_data[mat_id];
		    	
		    	// loop over the correct number of material cells
		    	for (int j=0; j<nmat_cells; ++j){
		    		these_material_data.push_back(correctly_typed_data[running_counter++]);
		    	}
		    	
		    }
		          
		    // add the material indices by keys
		    // the underlying api mat_add_celldata clears the existing data
		    for ( auto& kv: material_data){
		    	StateManager<MeshWrapper>::mat_add_celldata(field_name, kv.first, kv.second.data());
		    }
    	
    	} else if ( data_type == typeid(double)){
  		
  			// field data is doubles
  			
  			// allocate the data for unpacking
	      std::unordered_map<int,std::vector<double>> material_data;
	      
    	
				/////////////////////////////////////////////////////////
		    // We need to turn the flattened material cells into a correctly shaped
		    // ragged right structure for use as the material cells in the flat
		    // state wrapper. Just as in the flat state mesh field (and associated
		    // cell ids), we aren't removing duplicates, just concatnating by material
		    /////////////////////////////////////////////////////////
		    
		    
		    // reset the running counter
		    int running_counter=0;
		    
		    // loop over material ids on different nodes
		    for (int i=0; i<all_material_ids.size(); ++i){
		    
		    	// get the current working material
		    	int mat_id = all_material_ids[i];
		    	
		    	// get the current number of material cells for this material
		    	int nmat_cells = all_material_shapes[i];
		    	
		    	// get or create a reference to the correct material cell vector
		    	auto& these_material_data = material_data[mat_id];
		    	
		    	// loop over the correct number of material cells
		    	for (int j=0; j<nmat_cells; ++j){
		    		these_material_data.push_back(flat_data[running_counter++]);
		    	}
		    	
		    }
		          
		    // add the material indices by keys
		    // the underlying api mat_add_celldata clears the existing data
		    for ( auto& kv: material_data){
		    	StateManager<MeshWrapper>::mat_add_celldata(field_name, kv.first, kv.second.data());
		    }
    	
    	}
  	}
  }

  /*!
    @brief Get field stride
  */
  size_t get_field_stride(std::string field_name) {
  
  	// get the state vector base class shared pointer
    std::shared_ptr<StateVectorBase> pv =
        StateManager<MeshWrapper>::get(field_name);
         
    const std::type_info& data_type = pv->data_type(); 
    
    if ( data_type == typeid(double)) return 1;
    else if ( data_type == typeid(Wonton::Point<2>)) return 2;
    
  }

  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] index The index of the data field
    @return The Entity_kind enum for the entity type on which the field is
    defined
  */
  Entity_kind get_entity(std::string field_name) {
    return StateManager<MeshWrapper>::get(field_name)->get_kind();
  }
  
   /*!
    @brief Return the number of material cells for this node
    @return                    number of material cells for this node

    Return the number of material cells for this node. This number is the sum of
    the lengths of the material cell indices vectors. It is the number need to
    pass a flattened state vector.
  */ 
	int num_material_cells(){
		int n=0;
		for (auto& kv : StateManager<MeshWrapper>::material_cells_){
			n += kv.second.size();
		}
		return n;
	}
	
	/*!
    @brief Return the sorted vector of material ids actually used in the cell mat data.
    @return vector<int> of material id's (ordered)

    Return the sorted vector of material ids actually used in the cell mat data.
  */
  std::vector<int> const get_material_ids() const {

    std::vector<int> material_ids;
    for (auto& kv : StateManager<MeshWrapper>::material_cells_) material_ids.emplace_back(kv.first);
    
    std::sort(material_ids.begin(), material_ids.end());

    return material_ids;
  }

  /*!
    @brief Return the sorted vector of material shapes actually used in the cell mat data.
    @return vector<int> of material shapes (ordered)

    Return the sorted vector of material shapes actually used in the cell mat data.
  */
  std::vector<int> const get_material_shapes() const {

    std::vector<int> material_ids=get_material_ids();
    std::vector<int> material_shapes;
    
    for (auto m : material_ids) material_shapes.emplace_back(material_cells_.at(m).size());
    
    return material_shapes;
  }
  
  
  /*!
    @brief Return the sorted vector of material shapes actually used in the cell mat data.
    @return vector<int> of material shapes (ordered)

    Return the sorted vector of material shapes actually used in the cell mat data.
  */
  std::vector<int> const get_material_cells() const {

    std::vector<int> material_ids=get_material_ids();
    std::vector<int> material_cells;
    
    for (auto m : material_ids){
      std::vector<int> this_material_cells = StateManager<MeshWrapper>::get_material_cells(m);
      material_cells.insert(material_cells.end(), this_material_cells.begin(), this_material_cells.end());
    }
    
    return material_cells;
  }
	
};  // class Flat_State_Wrapper

}  // namespace Wonton

#endif  // WONTON_STATE_FLAT_FLAT_STATE_MM_WRAPPER_H_

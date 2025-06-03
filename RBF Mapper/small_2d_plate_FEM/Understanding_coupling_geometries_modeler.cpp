Understanding the coupling_geometries_modeler from the beginning

The mapper settings object is:
mMapperSettings : Parameters Object {
    "consistency_scaling": true,
    "destination_is_slave": true,
    "dual_mortar": false,
    "echo_level": 0,
    "linear_solver_settings": {},
    "modeler_name": "MappingGeometriesModeler",
    "modeler_parameters": {
        "destination_interface_sub_model_part_name": "Structure.DISPLACEMENT_left_edge",
        "destination_model_part_name": "destination",
        "is_interface_sub_model_parts_specified": true,
        "origin_interface_sub_model_part_name": "Structure.LineLoad2D_right_edge",
        "origin_model_part_name": "origin"
    },
    "precompute_mapping_matrix": false,
    "row_sum_tolerance": 1e-12,
    "solver_type": "skyline_lu_factorization"
}

/*
1) Everything starts with the MapperFactory class which is exposed to Python. From the Python level, you invoke the CreateMapper function with the following inputs:
- ModelPart& rModelPartOrigin 
- ModelPart& rModelPartDestination
- Parameters MapperSettings
The output type is: 
- Mapper<TSparseSpace, TDenseSpace>::Pointer
*/


static typename Mapper<TSparseSpace, TDenseSpace>::Pointer CreateMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters MapperSettings)
    {
        ModelPart& r_interface_model_part_origin = GetInterfaceModelPart(rModelPartOrigin, MapperSettings, "origin");
        ModelPart& r_interface_model_part_destination = GetInterfaceModelPart(rModelPartDestination, MapperSettings, "destination");

        KRATOS_ERROR_IF(!TSparseSpace::IsDistributed() && (r_interface_model_part_origin.IsDistributed() || r_interface_model_part_destination.IsDistributed())) << "Trying to construct a non-MPI Mapper with a distributed ModelPart. Please use \"CreateMPIMapper\" instead!" << std::endl;

        KRATOS_ERROR_IF(TSparseSpace::IsDistributed() && !r_interface_model_part_origin.IsDistributed() && !r_interface_model_part_destination.IsDistributed()) << "Trying to construct a MPI Mapper without a distributed ModelPart. Please use \"CreateMapper\" instead!" << std::endl;

        const std::string mapper_name = MapperSettings["mapper_type"].GetString();

        const auto& mapper_list = GetRegisteredMappersList();

        if (mapper_list.find(mapper_name) != mapper_list.end()) {
            // Removing Parameters that are not needed by the Mapper
            MapperSettings.RemoveValue("mapper_type");
            MapperSettings.RemoveValue("interface_submodel_part_origin");
            MapperSettings.RemoveValue("interface_submodel_part_destination");

            // TODO check why this works, Clone currently returns a unique ptr!!!
            return mapper_list.at(mapper_name)->Clone(r_interface_model_part_origin,
                                                      r_interface_model_part_destination,
                                                      MapperSettings);
        }
        else {
            std::stringstream err_msg;
            err_msg << "The requested Mapper \"" << mapper_name <<"\" is not not available!\n"
                    << "The following Mappers are available:" << std::endl;

            for (auto const& registered_mapper : mapper_list)
                err_msg << "\t" << registered_mapper.first << "\n";

            KRATOS_ERROR << err_msg.str() << std::endl;
        }
    }



    /* If the mapper_name specified in the MapperSettings object is implemented, it calls the Clone function associated
    with that specific mapper (in this case the Clone() method which belongs to the CouplingGeometryMapper)*/
    
    -------------------------------------------------------------------------------------------------------------------------
    The CreateMapper() is called from the class KratosMappingDataTransferOperator(CoSimulationDataTransferOperator) which is located in kratos_mapping.py, specifically from a function called
    ExecuteTransferData():
    
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        model_part_origin_name = from_solver_data.model_part_name
        variable_origin        = from_solver_data.variable
        identifier_origin      = from_solver_data.solver_name + "." + model_part_origin_name

        model_part_destination_name = to_solver_data.model_part_name
        variable_destination        = to_solver_data.variable
        identifier_destination      = to_solver_data.solver_name + "." + model_part_destination_name

        mapper_flags = self.__GetMapperFlags(transfer_options, from_solver_data, to_solver_data)

        identifier_tuple         = (identifier_origin, identifier_destination)
        inverse_identifier_tuple = (identifier_destination, identifier_origin)

        if identifier_tuple in self.__mappers:
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)
        elif inverse_identifier_tuple in self.__mappers:
            self.__mappers[inverse_identifier_tuple].InverseMap(variable_destination, variable_origin, mapper_flags)
        else:
            model_part_origin      = self.__GetModelPartFromInterfaceData(from_solver_data)
            model_part_destination = self.__GetModelPartFromInterfaceData(to_solver_data)

            if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
                mapper_create_fct = python_mapper_factory.CreateMPIMapper
            else:
                mapper_create_fct = python_mapper_factory.CreateMapper

            if self.echo_level > 0:
                info_msg  = "Creating Mapper:\n"
                info_msg += '    Origin: ModelPart "{}" of solver "{}"\n'.format(model_part_origin_name, from_solver_data.solver_name)
                info_msg += '    Destination: ModelPart "{}" of solver "{}"'.format(model_part_destination_name, to_solver_data.solver_name)

                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            mapper_creation_start_time = time()
            self.__mappers[identifier_tuple] = mapper_create_fct(model_part_origin, model_part_destination, self.settings["mapper_settings"].Clone()) # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other mappers

            if self.echo_level > 2:
                cs_tools.cs_print_info(colors.bold(self._ClassName()), "Creating Mapper took: {0:.{1}f} [s]".format(time()-mapper_creation_start_time,2))
            self.__mappers[identifier_tuple].Map(variable_origin, variable_destination, mapper_flags)
-------------------------------------------------------------------------------------------------------------------------------------------------------------
    


    /*2) The Clone() method belonging to the coupling_geometry_mapper is invoked (inside coupling_geometry_mapper.h).
    The input parameters are:
    - ModelPart& rModelPartOrigin 
    - ModelPart& rModelPartDestination
    - Parameters MapperSettings
    The output is an unique pointer to an instance of the coupling_geometry_mapper:
    - Kratos::make_unique<CouplingGeometryMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    Note: when the unique pointer is built, the constructor for coupling_geomtry_mapper is called with
    the parameters rModelPartOrigin, rModelPartDestination and JsonParameters
    */

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<CouplingGeometryMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    }

    /*3) The CouplingGeometryMapper constructor is called inside coupling_geometry_mapper.cpp:*/

    template<class TSparseSpace, class TDenseSpace>
CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CouplingGeometryMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters):
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mMapperSettings(JsonParameters)


        /* 4) An instance of the coupling_geometry modeler is created calling the method Create() in the 
        modeler_factory.cpp. The input is:
        - modeler_name
        - rModelPartOrigin
        - MapperSettings
        The output is:
        - A pointer to the instance of the coupling_geometry_modeler created
        */

        mpModeler = (ModelerFactory::Create(
        mMapperSettings["modeler_name"].GetString(),
        rModelPartOrigin.GetModel(),
        mMapperSettings["modeler_parameters"]));


        // The following piece of code is inside mapper_factory.cpp
        /// Checks if the modeler is registered
    typename Modeler::Pointer ModelerFactory::Create(
        const std::string& ModelerName, Model& rModel, const Parameters ModelParameters)
    {
        KRATOS_ERROR_IF_NOT(Has(ModelerName))
            << "Trying to construct a modeler: "
            << ModelerName << "\" which does not exist.\n"
            << "The available options (for currently loaded applications) are:\n"
            << KratosComponents< Modeler >() << std::endl;

        Modeler const& r_clone_modeler = KratosComponents< Modeler >::Get(ModelerName);
        return r_clone_modeler.Create(rModel, ModelParameters);
    }

    // Finally, the Create() method from the coupling geometry modeler is invoked
    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<MappingGeometriesModeler>(rModel, ModelParameters);
    }

    // Everywhere you see Pointer, remember that:

    #define KRATOS_CLASS_POINTER_DEFINITION(a) typedef Kratos::shared_ptr<a > Pointer; // This is a macro. If you say
    // KRATOS_CLASS_POINTER_DEFINITION(Modeler), it means typedef Kratos::shared_ptr<Modeler> Pointer

    /* 5) Adds the destination model part to the instance of the modeler created: */

    mpModeler->GenerateNodes(rModelPartDestination);

    // The method GenerateNodes() is implemented within mapping_geometries_modeler.h

    /// Adds the second model part to the modeler.
    void GenerateNodes(ModelPart& ThisModelPart) override
    {
        mpModels.push_back(&ThisModelPart.GetModel()); 
    }

    /* 6) The SetUpGeometryModel() method inside mapping_geometries_modeler.cpp is invoked: */

    mpModeler->SetupGeometryModel();

    // In the following steps, we will explain the logics behind the SetUpGeometryModel() method.

    /* 7) A new sub_model_part called coupling_model_part is created inside the origin model part */

    ModelPart& coupling_model_part = (mpModels[0]->HasModelPart("coupling"))
            ? mpModels[0]->GetModelPart("coupling")
            : mpModels[0]->CreateModelPart("coupling");

    /* 8) Create coupling conditions on the interface depending on the dimension, for both the origin and the destination domain */

    if (dim == 2)
        {
            CreateInterfaceLineCouplingConditions(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));

            //KRATOS_WATCH(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name))

            CreateInterfaceLineCouplingConditions(mpModels.back()->GetModelPart(destination_interface_sub_model_part_name));
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }

    /* 9) Understanding the CreateInterfaceLineCouplingConditions()*/

    /* 10) Creates a coupling_conditions model part inside the origin and destination interface */

     void MappingGeometriesModeler::CreateInterfaceLineCouplingConditions(ModelPart& rInterfaceModelPart)
    {
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_conditions = rInterfaceModelPart.GetSubModelPart("coupling_conditions");

    /* 11) It retrieves the parent model part (the interface model part is just a sub_model_part of this one) */
    
     const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();
     
     /* 12) Loop over the elements in the interface */
     
     for (size_t node_index = 0; node_index < rInterfaceModelPart.NumberOfNodes() - 1; ++node_index)
     
     /* 13) Retrives the id of the node in the interface and creates a vector of pointers to geometries*/
     
     interface_node_id = (rInterfaceModelPart.NodesBegin() + node_index)->Id();
     std::vector< GeometryPointerType> p_geom_vec;
     
     /* 14) It creates a vector of pointers to all the geometries which have nodes on the interface  (p_geom_vec) */
     
     for (auto& ele_it: root_mp.Elements())
            {
                auto p_geom = ele_it.pGetGeometry();
                for (size_t i = 0; i < p_geom->size(); i++)
                {
                    if ((*p_geom)[i].Id() == interface_node_id)
                    {
                        p_geom_vec.push_back(p_geom);
                    }
                }
            }
            
     Example for the origin model part:
     
     p_geom_vec.size() : 2
2 dimensional quadrilateral with four nodes in 2D space
    Working space dimension : 2
    Local space dimension   : 2

        Point 1  :  (2, 0, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 2  :  (2, 0.5, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 3  :  (1, 0.5, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 4  :  (1, 0, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Center   :  (1.5, 0.25, 0)


    Jacobian in the origin       : [2,2]((0,-0.5),(0.25,0))
2 dimensional quadrilateral with four nodes in 2D space
    Working space dimension : 2
    Local space dimension   : 2

        Point 1  :  (2, 0.5, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 2  :  (2, 1, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 3  :  (1, 1, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 4  :  (1, 0.5, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Center   :  (1.5, 0.75, 0)


    Jacobian in the origin       : [2,2]((0,-0.5),(0.25,0))
    
    /* 15) Loop over all geometries that have nodes on the interface */
    
    for (size_t interface_geom_index = 0; interface_geom_index < p_geom_vec.size(); interface_geom_index++)
            {
    GeometryType& r_interface_geom = *(p_geom_vec[interface_geom_index]);
    
   /* 16)  Loop over remaining interface nodes, see if any of them are nodes in the interface geom. If this is the case, create line coupling conditions between the points */
   
   for (size_t geom_node_index = 0; geom_node_index < r_interface_geom.size(); geom_node_index++)
                {
                    trial_geom_node_id = r_interface_geom[geom_node_index].Id();

                    for (size_t trial_index = node_index + 1; trial_index < rInterfaceModelPart.NumberOfNodes(); ++trial_index)
                    {
                        trial_interface_node_id = (rInterfaceModelPart.NodesBegin() + trial_index)->Id();
                        if (trial_geom_node_id == trial_interface_node_id)
                        {
                            // Another interface node was found in the same geom, make line condition between them
                            Geometry<GeometricalObject::NodeType>::PointsArrayType line_condition_points;
                            line_condition_points.push_back(rInterfaceModelPart.pGetNode(interface_node_id));
                            line_condition_points.push_back(rInterfaceModelPart.pGetNode(trial_interface_node_id));
                            coupling_conditions.CreateNewCondition("LineCondition2D2N", condition_id,
                                line_condition_points, rInterfaceModelPart.pGetProperties(0));
                            condition_id += 1;
                        }
                    }
                }
                
   Finally, you get:
   
   rInterfaceModelPart : -LineLoad2D_right_edge- model part
    Number of tables : 0
    Number of sub model parts : 2

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 3
        Number of Properties  : 1
        Number of Elements    : 0
        Number of Conditions  : 7
        Number of Constraints : 0

    -PLEASE_SPECIFY- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 3
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 3
            Number of Constraints : 0
    -coupling_conditions- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 0
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 2
            Number of Constraints : 0
    
    
   ----------------------------------------------------------------------------------------------------------------------------------------------------------
   rInterfaceModelPart : -DISPLACEMENT_left_edge- model part
    Number of tables : 0
    Number of sub model parts : 1

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 2
        Number of Properties  : 1
        Number of Elements    : 0
        Number of Conditions  : 1
        Number of Constraints : 0

    -coupling_conditions- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 0
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 1
            Number of Constraints : 0
            
    /* 17) Transfer everything into the coupling_model_part */
    
        ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        CopySubModelPart(coupling_interface_origin,
            mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));

        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");
        CopySubModelPart(coupling_interface_destination,
            mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));
            
        Finally, coupling_model_part is:
        
        coupling_model_part : -coupling- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 2
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0

    -interface_destination- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 2
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 1
            Number of Constraints : 0
    -interface_origin- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 3
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 2
            Number of Constraints : 0
            
     /* 18) If working_dim==2 & local_dim==1, find the intersection between the nodes in both interfaces */
     
        if (working_dim == 2 && local_dim == 1)
        {
            MappingIntersectionUtilities::FindIntersection1DGeometries2D(
                coupling_interface_origin,
                coupling_interface_destination,
                coupling_model_part, 1e-6);
 
    
     Now, we will explain this MappingIntersectionUtilities::FindIntersection1DGeometries2D
     
     /*  19) Iterate over all the combinations of conditions in the coupling_interface_origin and coupling_interface_destination submodel parts */
     
     for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
        condition_a_itr != rModelPartDomainA.ConditionsEnd();
        ++condition_a_itr)
    {
        for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
            condition_b_itr != rModelPartDomainB.ConditionsEnd();
            ++condition_b_itr)
        {
            if (FindOverlapExtents1DGeometries2D(condition_a_itr->GetGeometry(), condition_b_itr->GetGeometry(), dummy))
            {
                rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                    condition_a_itr->pGetGeometry(), condition_b_itr->pGetGeometry()));
            }
        }
    }
    
    - If an overlap is found between the 2 geometries (return==true), it adds a CouplingGeometry to the coupling_model_part (note that the CouplingGeometry constructor is called with  condition_a_itr->pGetGeometry() and condition_b_itr->pGetGeometry() )
    
     /// Constructor for coupling one master to one slave geometry.
    CouplingGeometry(
        GeometryPointer pMasterGeometry,
        GeometryPointer pSlaveGeometry)
        : BaseType(PointsArrayType(), &(pMasterGeometry->GetGeometryData()))
    {
        mpGeometries.resize(2);

        mpGeometries[0] = pMasterGeometry;
        mpGeometries[1] = pSlaveGeometry;
    }
    
    ----------------------------------------------------------------------------------------------
    
    The final output for the coupling_model_part is:
    coupling_model_part : -coupling- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 2
    Current solution step index : 0

    Number of Geometries  : 2
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0

    -interface_destination- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 2
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 1
            Number of Constraints : 0
    -interface_origin- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 3
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 2
            Number of Constraints : 0
            
   -------------------------------------------------------------------------------------------------
   If you print the coupling geometries, you get:
   
   geometry_itr : Coupling geometry that holds a master and a set of slave geometries.
    Working space dimension : 2
    Local space dimension   : 1

        Center   :  (2, 0.25, 0)


    CouplingGeometry with 2 geometries.
r_geom_master : 1 dimensional line in 2D space
    Working space dimension : 2
    Local space dimension   : 1

        Point 1  :  (2, 0.5, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 2  :  (2, 0, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Center   :  (2, 0.25, 0)


    Jacobian     : [2,1]((0),(-0.25))
r_geom_slave : 1 dimensional line in 2D space
    Working space dimension : 2
    Local space dimension   : 1

        Point 1  :  (2, 1, 0)
    Dofs :
        Fix DISPLACEMENT_X degree of freedom
        Fix DISPLACEMENT_Y degree of freedom
        Fix DISPLACEMENT_Z degree of freedom

        Point 2  :  (2, 0, 0)
    Dofs :
        Fix DISPLACEMENT_X degree of freedom
        Fix DISPLACEMENT_Y degree of freedom
        Fix DISPLACEMENT_Z degree of freedom

        Center   :  (2, 0.5, 0)


    Jacobian     : [2,1]((0),(-0.5))
    
    
    
    
    geometry_itr : Coupling geometry that holds a master and a set of slave geometries.
    Working space dimension : 2
    Local space dimension   : 1

        Center   :  (2, 0.25, 0)


    CouplingGeometry with 2 geometries.
r_geom_master : 1 dimensional line in 2D space
    Working space dimension : 2
    Local space dimension   : 1

        Point 1  :  (2, 0.5, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Point 2  :  (2, 0, 0)
    Dofs :
        Free DISPLACEMENT_X degree of freedom
        Free DISPLACEMENT_Y degree of freedom
        Free DISPLACEMENT_Z degree of freedom

        Center   :  (2, 0.25, 0)


    Jacobian     : [2,1]((0),(-0.25))
r_geom_slave : 1 dimensional line in 2D space
    Working space dimension : 2
    Local space dimension   : 1

        Point 1  :  (2, 1, 0)
    Dofs :
        Fix DISPLACEMENT_X degree of freedom
        Fix DISPLACEMENT_Y degree of freedom
        Fix DISPLACEMENT_Z degree of freedom

        Point 2  :  (2, 0, 0)
    Dofs :
        Fix DISPLACEMENT_X degree of freedom
        Fix DISPLACEMENT_Y degree of freedom
        Fix DISPLACEMENT_Z degree of freedom

        Center   :  (2, 0.5, 0)


    Jacobian     : [2,1]((0),(-0.5))
    
    /* 20) Creates the quadrature points */
    
    MappingIntersectionUtilities::CreateQuadraturePointsCoupling1DGeometries2D(
                coupling_model_part, 1e-6);
                
    /* 21) Iterate over the coupling geometries inside couplig_model_part and for each coupling geometry, creates a vector with the overlap extents */
    
    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
         auto& r_geom_master = geometry_itr->GetGeometryPart(0);
        auto& r_geom_slave = geometry_itr->GetGeometryPart(1);
        CoordinatesArrayType local_parameter_1 = ZeroVector(3);
        CoordinatesArrayType local_parameter_2 = ZeroVector(3);

        std::vector<array_1d<double, 3>> overlap_extents;
        
    /* 22) Find the overlap extents between the master and slave in order to set the integration points */
    
    KRATOS_ERROR_IF_NOT(MappingIntersectionUtilities::FindOverlapExtents1DGeometries2D(
            r_geom_master, r_geom_slave, overlap_extents, 1e-6))
            << "Lines do not intersect." << std::endl;
            
     The output is overlap_extents in the geometrical space:
     overlap_extents : [[3](2,0.5,0), [3](2,0,0)]
    
    /* 23) Find the local coordinates (natural coordinates) of the overlap extents */
    
    r_geom_master.PointLocalCoordinates(local_parameter_1, overlap_extents[0]); // min of overlap
    r_geom_master.PointLocalCoordinates(local_parameter_2, overlap_extents[1]); // max of overlap
    
    The output for the last case is: 
    local_parameter_1 : [3](-1,0,0)
    local_parameter_2 : [3](1,0,0)
    
    /* 24) For every coupling geomtry, creates a vector of integration points anc call the IntegrationPoints1D method*/
    
    IntegrationPointsArrayType integration_points(IntegrationPointsPerSpan);
    
    typename IntegrationPointsArrayType::iterator integration_point_iterator = integration_points.begin();

    IntegrationPointUtilities::IntegrationPoints1D(
            integration_point_iterator,
            IntegrationPointsPerSpan,
            local_parameter_1[0], local_parameter_2[0]);
    
    The output for the first coupling geometry is: 
    integration_points : [3 dimensional integration point(-0.57735 , 0 , 0), weight = 1, 3 dimensional integration point(0.57735 , 0 , 0), weight = 1]
    
    /* 25) Determine quadrature point locations of span on master and then create quadrature point geometries */
    
        GeometriesArrayType quadrature_point_geometries_master(IntegrationPointsPerSpan);
        CreateQuadraturePointsUtility<NodeType>::Create(
            r_geom_master, quadrature_point_geometries_master, integration_points, 1); // 1 is the number of shape functions derivatives 
        
    	The Create() method inside quadrature_points_utility.h: 
    	
    	/// creates a quadrature point geometry on a provided location.
        static void Create(
            GeometryType& rGeometry,
            typename GeometryType::GeometriesArrayType& rResultGeometries,
            typename GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
            SizeType NumberOfShapeFunctionDerivatives)
        {
            KRATOS_ERROR_IF(NumberOfShapeFunctionDerivatives > 1)
                << "Create can only compute shape functions up to an derivative order of 1. "
                << "Demanded derivative order: " << NumberOfShapeFunctionDerivatives << std::endl;

            // Resize containers.
            if (rResultGeometries.size() != rIntegrationPoints.size())
                rResultGeometries.resize(rIntegrationPoints.size());

            auto default_method = rGeometry.GetDefaultIntegrationMethod();

            Vector N;
            Matrix DN_De;
            for (IndexType i = 0; i < rIntegrationPoints.size(); ++i)
            {
                rGeometry.ShapeFunctionsValues(N, rIntegrationPoints[i]);

                Matrix N_matrix = ZeroMatrix(1, N.size());
                if (NumberOfShapeFunctionDerivatives >= 0) {
                    for (IndexType j = 0; j < N.size(); ++j)
                    {
                        N_matrix(0, j) = N[j];
                    }
                }

                /// Get Shape Function Derivatives DN_De, ...
                if (NumberOfShapeFunctionDerivatives > 0) {
                    rGeometry.ShapeFunctionsLocalGradients(DN_De, rIntegrationPoints[i]);
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    default_method, rIntegrationPoints[i],
                    N_matrix, DN_De);

                rResultGeometries(i) = CreateQuadraturePointsUtility<TPointType>::CreateQuadraturePoint(
                    rGeometry.WorkingSpaceDimension(), rGeometry.LocalSpaceDimension(),
                    data_container, rGeometry);
            }
        }
        
        - For each integration point, it creates a QuadraturePoint which contains all the information of the shape functions and it derivatives at the integration point 
        
        
        /* 26) Transfer the quadrature points to slave (projection step) */
        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i)
        {
            CoordinatesArrayType local_parameter_slave = ZeroVector(3);
            r_geom_slave.PointLocalCoordinates(local_parameter_slave, quadrature_point_geometries_master[i].Center());

            KRATOS_WATCH(local_parameter_slave)

            integration_points[i].X() = local_parameter_slave[0];
        }
        
        Nots:
        - The quadrature_point_geometries_master[i].Center() method gives the position of the i_th quadrature point geometry in the physical space of the master domain 
        - r_geom_slave.PointLocalCoordinates(local_parameter_slave, quadrature_point_geometries_master[i].Center()) gives the natural coordinates of the quadrature points in the slave parameter space 
        
        For example, for the first two quadrature points geometries:
        quadrature_point_geometries_master[i].Center() : Point (2, 0.394338, 0) // Coordinate of the QP in the physical space 
local_parameter_slave : [3](0.211325,0,0) // Natural coordinate of the QP
quadrature_point_geometries_master[i].Center() : Point (2, 0.105662, 0)
local_parameter_slave : [3](0.788675,0,0)

integration_points : [3 dimensional integration point(0.211325 , 0 , 0), weight = 1, 3 dimensional integration point(0.788675 , 0 , 0), weight = 1]

	quadrature_point_geometries_master[i].Center() : Point (2, 0.894338, 0)
local_parameter_slave : [3](-0.788675,0,0)
quadrature_point_geometries_master[i].Center() : Point (2, 0.605662, 0)
local_parameter_slave : [3](-0.211325,0,0)

integration_points : [3 dimensional integration point(-0.788675 , 0 , 0), weight = 1, 3 dimensional integration point(-0.211325 , 0 , 0), weight = 1]

	/* 27) Create slave quadrature point geometries */
	emplate<class TSparseSpace, class TDenseSpace>
CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CouplingGeometryMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters):
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mMapperSettings(JsonParameters)
{
    JsonParameters.ValidateAndAssignDefaults(GetMapperDefaultSettings(
	GeometriesArrayType quadrature_point_geometries_slave(IntegrationPointsPerSpan);
        CreateQuadraturePointsUtility<NodeType>::Create(
            r_geom_slave, quadrature_point_geometries_slave, integration_points, 1);
            
       
       /* 28) Add the Quadrature point geometry conditions to the result model part */
       
       const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
            ? 1
            : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
        for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i) {
            rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
        }
        
        The output is:
        
        rModelPartCoupling : -coupling- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 2
    Current solution step index : 0

    Number of Geometries  : 2
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 4  // this is the important change: 4 quadrature point geometry conditions!
        Number of Constraints : 0

    -interface_destination- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 2
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 1
            Number of Constraints : 0
    -interface_origin- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 3
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 2
            Number of Constraints : 0
            
            
            
       // Now we come back to the coupling_geometry_mapper.cpp file
       
       /* 29) It constructs two interface vector containers */
       
       mpInterfaceVectorContainerMaster = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceMaster);
    	mpInterfaceVectorContainerSlave = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceSlave);
    	
    	/* 30) It creates the linear solver */
    	
    	
    	    this->CreateLinearSolver();
    	    
    	    Here is the CreateLinearSolver() implementation:
    	
    	
    	template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    if (mMapperSettings["linear_solver_settings"].Has("solver_type")) {
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mMapperSettings["linear_solver_settings"]);
    }
    else {
        // TODO - replicate 'get fastest solver'
        mMapperSettings.AddString("solver_type", "skyline_lu_factorization");
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mMapperSettings);
    }
}

	/* 31) Call the initialize interface method */
	
	
	 this->InitializeInterface();
	 
	 template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{

	/* 32) Creates 2 instances of the class CouplingGeometryLocalSystem */
	
	CouplingGeometryLocalSystem ref_projector_local_system(nullptr, true, dual_mortar, direct_map_to_destination);
    	CouplingGeometryLocalSystem ref_slave_local_system(nullptr, false, dual_mortar, direct_map_to_destination);
        
        /* This class can be found in coupling_geometry_mapper.h */
        class CouplingGeometryLocalSystem : public MapperLocalSystem
{
public:

    explicit CouplingGeometryLocalSystem(GeometryPointerType pGeom,
                                         const bool IsProjection,
                                         const bool IsDualMortar,
                                         const bool IsDestinationIsSlave
                                         )
        : mpGeom(pGeom),
          mIsProjection(IsProjection),
          mIsDualMortar(IsDualMortar),
          mIsDestinationIsSlave(IsDestinationIsSlave)
        {}
        
        void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;
                      
                      
        /* 33) Creates the local systems objects which will finally give rise to the mapping matrices (Note: 1 local matrix per gaussian point) */
        
        MapperUtilities::CreateMapperLocalSystemsFromGeometries(ref_projector_local_system,
                             mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystemsProjector);

    	MapperUtilities::CreateMapperLocalSystemsFromGeometries(ref_slave_local_system,
                             mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystemsSlave);

        
        /* 34) Assign interface equations id´s */
        
        AssignInterfaceEquationIds(); // Has to be done every time in case of overlapping interfaces!
        
        /* 35) Creates the projector interface mass matrix - interface_matrix_projector */
        
        const std::size_t num_nodes_interface_slave = mpCouplingInterfaceSlave->NumberOfNodes();
    	const std::size_t num_nodes_interface_master = mpCouplingInterfaceMaster->NumberOfNodes();
    	mpMappingMatrix = Kratos::make_unique<MappingMatrixType>(num_nodes_interface_slave, num_nodes_interface_master);
    	
    	// If you dereference and print mpMappingMatrix, you get:
    	
    	*mpMappingMatrix : [2,3]((0,0,0),(0,0,0))
    	
    	/* 36) It builds the slave mapping matrix (int_gammaD N_D * N_D dgamma_D) */
    	
    	MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrix(
        mpMappingMatrixSlave,
        mpInterfaceVectorContainerSlave->pGetVector(),
        mpInterfaceVectorContainerSlave->pGetVector(),
        mpInterfaceVectorContainerSlave->GetModelPart(),
        mpInterfaceVectorContainerSlave->GetModelPart(),
        mMapperLocalSystemsSlave,
        0); // The echo-level is no longer neeed here, refactor in separate PR
        
        Note: Check the BuildMappingMatrix() method inside mapping_matrix_utilities.cpp
        
        The output is: 
        
        *mpMappingMatrixSlave : [2,2]((0.666667,0.333333),(0.333333,0.666667))
        
        
        /* 37) It builds the projector mapping matrix (int_gammaDO N_D * N_O dgamma_DO) */
        
        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrix(
        mpMappingMatrixProjector,
        mpInterfaceVectorContainerMaster->pGetVector(),
        mpInterfaceVectorContainerSlave->pGetVector(),
        mpInterfaceVectorContainerMaster->GetModelPart(),
        mpInterfaceVectorContainerSlave->GetModelPart(),
        mMapperLocalSystemsProjector,
        0); // The echo-level is no longer neeed here, refactor in separate PR
        
        The output is:
        *mpMappingMatrixProjector : [2,3]((0.416667,0.5,0.0833333),(0.0833333,0.5,0.416667))
        
        /* 38) Then, if requested, it performs consistency scaling */
        
        if (mMapperSettings["consistency_scaling"].GetBool()) {
        EnforceConsistencyWithScaling(*mpMappingMatrixSlave, *mpMappingMatrixProjector, 1.1);
    	}
        
        
    	
        

     
     
    

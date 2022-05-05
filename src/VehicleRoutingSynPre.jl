module VehicleRoutingSynPre

using JLD2, Random
# Write your package code here.
include("ParticleSwarm.jl")

export Particle,
       load_data,
       generate_empty_particle,
       generate_particles,
       find_group_of_node,
       find_starttime,
       initial_insert_service,
       find_service_request,
       find_SYN,
       find_PRE,
       find_vehicle_service,
       feasibility,
       compatibility,
       find_compat_vehicle_node,
       create_empty_slot,
       find_vehicle_by_service,
       example,
       insert_vehicle_to_service,
       check_assigned_node,
       find_remain_service,
       total_distance,
       tardiness,
       objective_value,
       insert_PRE,
       complete,
       find_location_by_node_service,
       generate_example,
       example2
end

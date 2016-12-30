bin_PROGRAMS = netjob
INCLUDES = -I/usr/local/include/igraph

netjob_LDADD= -ligraph -lgsl -lgslcblas
#-lboost_system -lboost_thread
netjob_LDFLAGS= -static
netjob_SOURCES=	main.cc\
				constant.h\
				neighbor_info.h\
				network.cc\
				network.h\
				parameter.cc\
				parameter.h\
				pearl_network.h\
				runner.h\
				runner_manager.h\
				static_network.cc\
				static_network.h\
				dynamic_network.cc\
				dynamic_network.h\
				pearl_network.cc\
				pearl_network.h\
				step_test.cc\
				step_test.h\
				temp_node_data.h\
				undirected_graph.cc\
				undirected_graph.h\
				util_boost.cc\
				util_boost.h\
				util.cc\
				util_constant.cc\
				util_constant.h\
				util.h\
				util_igraph.cc\
				util_igraph.h\
				util_gsl.cc\
				util_gsl.h\
				disease_dynamics.cc\
				disease_dynamics.h\
				simulator.h\
				simulator.cc\
				source_identification.h\
				source_identification.cc\
				bfs_tree.h\
				bfs_tree.cc\
				step_source_identification.h\
				step_source_identification.cc\
				step_source_identification_v2.h\
				step_source_identification_v2.cc\
				step_get_small_component.h\
				step_get_small_component.cc\
				step_calc_static_feature.h\
				step_calc_static_feature.cc\
				step_evolve_large_component.h\
				step_evolve_large_component.cc\
				step_reindex_gml.h\
				step_reindex_gml.cc\
				step_select_nelder_mead_solutions.h\
				step_select_nelder_mead_solutions.cc\
				step_simulate_evolve.h\
				step_simulate_evolve.cc\
				step_determine_relationship.h\
				step_determine_relationship.cc\
				step_gen_contact_parts.h\
				step_gen_contact_parts.cc

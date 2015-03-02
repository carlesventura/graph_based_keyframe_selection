/*
 * graph_properties.hpp
 *
 *  Created on: Feb 28, 2012
 *      Author: guillem
 */

#ifndef GRAPH_PROPERTIES_HPP_
#define GRAPH_PROPERTIES_HPP_

#include <imageplus/core.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace imageplus {
	namespace math {
		namespace graphs {

		/*
		 * Graph types
		 */

		//! Default type for the graphs. Bidirectional (two edges to represent an undirected connection) and possibility to acces to the node incoming edges
		typedef boost::bidirectionalS 	kGraphBidirectional; // better than directed, but uses twice the space

		//! Undirected type
		typedef boost::undirectedS		kGraphUndirected;

		//! Directed type. cannot acces the incoming edge connections
		typedef boost::directedS		kGraphDirected;

		/*
		 * Redefine it if you change the graph
		 */
		//! Default graph type
		typedef kGraphBidirectional DefaultGraphType;

		/*********************************
		 *         Default Property Maps *
		 *********************************/

		//! Auxiliary type
		typedef boost::adjacency_list_traits < boost::listS, boost::vecS, DefaultGraphType > Traits;

		/*
		 * Default node properties
		 */

		//! Default structure for node properties
		struct NodeDefaultProperties {
			//! Node name
			std::string name;

			//! Node id
			uint64 id;

			//! Default constructor
			NodeDefaultProperties() : name(""), id(0) {}
		};

		/*
		 * Properties used by the max flow algorithm. mandatory.
		 */

		//! Properties for the maxflow algorithms
		struct NodeMaxFlowProperties : public NodeDefaultProperties {
			//! auxiliary types
			boost::default_color_type color;

			//! auxiliary types
			float64 distance;

			//! auxiliary types
			Traits::edge_descriptor	predecessor;

			//! Default constructor
			NodeMaxFlowProperties() : NodeDefaultProperties(), color(boost::white_color), predecessor() {}
		};

		/*
		 * Default edge properties
		 */
		//! Default Edge properties
		struct EdgeDefaultProperties {

			//! Weight of the edge
			float64 weight;

			//! Id of the edge
			float64 id;

			//! Default constructor
			EdgeDefaultProperties() : weight(0), id(0) {}

			//! Check if two edges are equal. used in the BoostGraph class
			bool operator==(const EdgeDefaultProperties& p) {
				return (p.weight == weight && id == p.id);
			}
		};

		/*
		 * Properties
		 */

		//! Maxflow edge properties
		struct EdgeMaxFlowProperties : public EdgeDefaultProperties {

			//! Capacity of the edge
			float64 capacity;

			//! Residual capacity
			float64 residual_capacity;

			//! auxiliary type. MaxFlowMinCut class handles this automatically
			Traits::edge_descriptor reverse_edge;

			//! Default constructor
			EdgeMaxFlowProperties() : EdgeDefaultProperties(), capacity(0), residual_capacity(0), reverse_edge() {}

			//! Check if two edges are equal. used in the BoostGraph class
			bool operator==(const EdgeMaxFlowProperties& p) {
				return (EdgeDefaultProperties::operator==(p) && capacity == p.capacity && residual_capacity == p.residual_capacity);
			}
		};

		}
	}
}


#endif /* GRAPH_PROPERTIES_HPP_ */

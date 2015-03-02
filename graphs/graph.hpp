// --------------------------------------------------------------
// Copyright (C)
// Universitat Politecnica de Catalunya (UPC) - Barcelona - Spain
// --------------------------------------------------------------

/*!
 * \file    graph.hpp
 * \author  Guillem Palou
 * \date    Oct 20, 2011
 */

#ifndef IMAGEPLUS_MATH_GRAPHS_GRAPH_HPP
#define IMAGEPLUS_MATH_GRAPHS_GRAPH_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
//#include <boost/graph/subgraph.hpp>
#include <imageplus/math/graphs/graph_properties.hpp>

namespace imageplus
{
namespace math
{
namespace graphs
{
    //! \brief Wrapper class for boost graph
    //! \warning Do not change the default properties unless you know what you are doing
    template<class DirectedS = DefaultGraphType, class NodeProperties = NodeDefaultProperties, class EdgeProperties = EdgeDefaultProperties> // should be boost::directedS, boost::undirectedS, or boost::bidirectionalS (default)
    class BoostGraph
    {
    public:
        //! an adjacency_list boost graph. basic structure
        typedef typename 	boost::adjacency_list<
                            boost::listS, // edge container allow parallel edges
                            boost::vecS, // vertex container
                            DirectedS,
                            NodeProperties,
                            EdgeProperties
                            > GraphContainer;

        //! Node Properties
        typedef NodeProperties	NodePropertiesType;

        //!Edge properties
        typedef EdgeProperties  EdgePropertiesType;

        // There are only the most typical ones

        //! Node type
        typedef typename boost::graph_traits<GraphContainer>::vertex_descriptor 				Node;

//		typedef	typename boost::subgraph< boost::adjacency_list< boost::listS, boost::vecS, DirectedS, NodeProperties, EdgeProperties > > SubGraph;


        //! Edge type
        typedef typename boost::graph_traits<GraphContainer>::edge_descriptor 					Edge;

        //! Auxiliary types
        typedef typename std::pair<Node, Node> 													NodePair;

        //! Auxiliary types
        typedef typename std::pair<Edge, Edge> 													EdgePair;

        // Iterators

        //! boost iterator for all nodes of the graph
        typedef typename boost::graph_traits<GraphContainer>::vertex_iterator 		boost_node_iterator;

        //! boost iterator for adjacency nodes
        typedef typename boost::graph_traits<GraphContainer>::adjacency_iterator 	boost_adjacency_iterator;

        //! boost iterator for out edges of a node
        typedef typename boost::graph_traits<GraphContainer>::out_edge_iterator 	boost_out_edge_iterator;

        //! boost iterator for in edges of a node
        typedef typename boost::graph_traits<GraphContainer>::in_edge_iterator 	boost_in_edge_iterator;

        // Properties

        // constructors etc.
        //! default constructor
        BoostGraph()
        {}

        /*!
         * Constructor taking the number of nodes
         *
         * \param[in] N : number of graph nodes;
         */
        BoostGraph(uint32 N)
        :   _graph(N)
        {}

        /*!
         * Copy constructor
         * \param[in] g : another BoostGraph
         */
        BoostGraph(const BoostGraph& g)
        :   _graph(g._graph)
        {}

        virtual
        ~BoostGraph()
        {}

        // Node Information methods

        //! \return The number of the nodes in the graph
        uint64 num_nodes()
        {
            return boost::num_vertices(_graph);
        }

        // Edge Information methods

        //! \return The number of the edges in the graph
        uint64 num_edges()
        {
            return boost::num_edges(_graph);
        }

        //! \return The number of the edges going out of a node
        //! \param[in] n : The Node to get the going out edges
        uint64 num_out_edges(Node n)
        {
            return boost::out_degree(n, _graph);
        }

        //! \return The number of edges going inside (entering) a node
        //! \param[in] n : The Node to get the entering edges
        uint64 num_in_edges(Node n)
        {
            return boost::in_degree(n, _graph);
        }

        //! \return The source A of an edge A->B
        //! \param e : The Edge A->B
        Node source(Edge e)
        {
            return boost::source(e, _graph);
        }

        //! \return The target B of an edge A->B
        //! \param e : The Edge A->B
        Node target(Edge e)
        {
            return boost::target(e, _graph);
        }

        //! \return true if exists an edge between Node A and Node B
        //! \param[in] a : node A (A->B)
        //! \param[in] b : node B (A->B)
        bool edge_exists(Node a, Node b)
        {
            return boost::edge(a,b,_graph).second;
        }

        //! \return Edge from A to B
        //! \param[in] a : node A (A->B)
        //! \param[in] b : node B (A->B)
        Edge edge(Node a, Node b)
        {
            return boost::edge(a,b,_graph).first;
        }

        /*
         * structure modification methods
         */

        //! \brief Erases all nodes and edges
        void clear()
        {
            _graph.clear();
        }

        //! \brief Adds a node to the graph
        //! \returns Reference to the node
        Node add_node()
        {
            return boost::add_vertex(_graph);
        }

        //! \brief Removes a node A from the graph
        void remove_node(Node a)
        {
            boost::clear_vertex(a,_graph);
            boost::remove_vertex(a,_graph);
        }

        //! Adds an edge to the graph
        //! \param[in] e : reference to an edge. use with caution
        //! \warning Be careful with this function. The edge to be added is an edge of the same graph (deleted a priori, for example)
        //! \returns the added edge
        Edge add_edge(Edge e)
        {
            return boost::add_edge(boost::source(e,_graph), boost::target(e,_graph), edge_properties(e), _graph).first;
        }

        //! Adds an edge to the graph
        //! \param[in] a : source node
        //! \param[in] b : target node
        //! \param[in] p : Properties of the edge
        //! \returns the added edge
        Edge add_edge(Node a, Node b, EdgePropertiesType p)
        {
            return boost::add_edge(a, b, p, _graph).first;
        }

        //! Adds an edge to the graph
        //! \param[in] a : source node
        //! \param[in] b : target node
        //! \returns the added edge
        Edge add_edge(Node a, Node b)
        {
            return boost::add_edge(a, b, EdgePropertiesType(), _graph).first;
        }

        //! Adds an edge to the graph
        //! \param[in] e : reference to an edge. use with caution
        void remove_edge(Edge e)
        {
            boost::remove_edge(e,_graph);
        }

        //! Removes an edge of the graph
        //! \param[in] a : source node
        //! \param[in] b : target node
        void remove_edge(Node a, Node b)
        {
            boost::remove_edge(a,b,_graph);
        }

        //! Removes an edge with given properties
        //! \param[in] a : source node
        //! \param[in] b : target node
        //! \param[in] p : properties
        void remove_edge(Node a, Node b, EdgePropertiesType& p)
        {
            for (out_edge_iterator e = out_edges_begin(a); e != out_edges_end(a); ++e) {
                if (edge_properties(*e) == p && target(*e) == b) {
                    remove_edge(*e);
                    break;
                }
            }
        }

        //! Removes an edge and its reverse of the graph
        //! \param[in] a : source node
        //! \param[in] b : target node
        void remove_undirected_edge(Node a, Node b)
        {
            boost::remove_edge(a,b,_graph);
            boost::remove_edge(b,a,_graph);
        }

        //! \return the boost graph structure (used for the algorithms) (adjacency list)
        GraphContainer& graph()
        {
            return _graph;
        }

        //! \return the boost graph structure const (used for the algorithms) (adjacency list)
        const GraphContainer& graph() const
        {
            return _graph;
        }

        /*
         * Properties modifiers
         */

        //! \param[in] n : node to get the Properties from
        //! \return NodeProperties structure for the node n
        NodeProperties& node_properties(Node n)
        {
            return _graph[n];
        }

        //! \param[in] e : the Edge to get the propierties from
        //! \return EdgeProperties structure for the edge e
        EdgeProperties& edge_properties(Edge e)
        {
            return _graph[e];
        }

        /*!
         * \returns the connected components of the graph.
         *
         * See http://www.boost.org/doc/libs/1_55_0/libs/graph/doc/connected_components.html
         */
		std::vector<int> connected_components()
		{
			std::vector<int> component(boost::num_vertices(_graph));
			boost::connected_components(_graph, &component[0]);
			return component;
		}

        /*!
         * Wrapper to the boost iterators
         * Making them look more like the STL iterators
         */

        /*!
         * \cond SKIP_DOC
         *
         * Private common implementation. It works as an ordinary STL iterator.  The only thing you should use is the de-reference operator
         */
        template<class boost_iterator_type>
        class node_iterator_base : public std::iterator<std::forward_iterator_tag, Node>
        {
            //! auxiliary typedef
            typedef typename std::pair<boost_iterator_type,boost_iterator_type>			node_iterator_pair;

        public:

            //!
            node_iterator_base(node_iterator_pair p) : it(p) {
            }

            node_iterator_base(boost_iterator_type a, boost_iterator_type b) : it(a,b) {
            }

            node_iterator_base& operator++() {
                ++it.first;
                return *this;
            }

            bool operator==(node_iterator_base& v) {
                return *it.first == *v.it.first;
            }

            bool operator!=(const node_iterator_base& v) const {
                return *it.first != *v.it.first;
            }

            //! Dereference node iterator
            //! \return
            Node operator*() {
                return (*it.first);
            }

        private:

            node_iterator_pair it;
        };
        /*!
         * \endcond SKIP_DOC
         */

        /*Types of node iterators*/

        //! Iterator of all the nodes in the graph
        typedef node_iterator_base<boost_node_iterator> 		node_iterator;

        //! Iterator of adjacent nodes of a given node
        typedef node_iterator_base<boost_adjacency_iterator> 	adjacent_node_iterator;

        //! Beginning of the nodes in the graph
        //! \return node iterator pointing to the first node in the graph
        node_iterator nodes_begin()
        {
            return node_iterator(boost::vertices(_graph));
        }

        //! ending of the nodes in the graph
        //! \return node iterator pointing to the end of the graph
        node_iterator nodes_end()
        {
            return node_iterator(boost::vertices(_graph).second, boost::vertices(_graph).second);
        }

        //! Beginning of adjacent nodes of a given node in the graph
        //! \param[in] n : node to search its adjacent nodes for
        //! \return node iterator pointing to the first node in the graph
        adjacent_node_iterator adjacent_nodes_begin(Node n)
        {
            return adjacent_node_iterator(boost::adjacent_vertices(n,_graph));
        }

        //! Beginning of adjacent nodes of a given node in the graph
        //! \param[in] n : node to search its adjacent nodes for
        //! \return node iterator pointing to the first node in the graph
        adjacent_node_iterator adjacent_nodes_end(Node n)
        {
            return adjacent_node_iterator(boost::adjacent_vertices(n,_graph).second, boost::adjacent_vertices(n,_graph).second);
        }

        /*!
         * Wrapper to the Edge iterator class
         */

        /*!
         * \cond SKIP_DOC
         *
         * Private common implementation. It works as an ordinary STL iterator. The only thing you should use is the de-reference operator
         */
        template<class boost_iterator_type>
        class edge_iterator_base : public std::iterator<std::forward_iterator_tag, Edge>
        {
            typedef typename std::pair<boost_iterator_type,boost_iterator_type>					edge_iterator_pair;

            /*typedef typename boost::graph_traits<GraphContainer>::out_edge_iterator 			boost_out_edge_iterator;
            typedef typename std::pair<boost_out_edge_iterator,boost_out_edge_iterator>			out_edge_iterator_pair;*/

        public:
            edge_iterator_base(edge_iterator_pair e) : it(e) {}
            edge_iterator_base(boost_iterator_type a,boost_iterator_type b) : it(a,b) {}

            edge_iterator_base& operator++() {
                ++it.first;
                return *this;
            }

            bool operator==(edge_iterator_base& v) {
                return *it.first == *v.it.first;
            }

            bool operator!=(const edge_iterator_base& v) const {
                return *it.first != *v.it.first;
            }

            Edge operator*() {
                return (*it.first);
            }
        private:
            edge_iterator_pair it;
        };
        /*!
         * \endcond SKIP_DOC
         */


        //! Iterator for all the edges in the graph
        typedef edge_iterator_base<typename boost::graph_traits<GraphContainer>::edge_iterator> 		edge_iterator;

        //! Iterator for all the out edges of a given node
        typedef edge_iterator_base<typename boost::graph_traits<GraphContainer>::out_edge_iterator> 	out_edge_iterator;

        //! Iterator for all the in edges of a given node
        typedef edge_iterator_base<typename boost::graph_traits<GraphContainer>::in_edge_iterator> 		in_edge_iterator;

        /*Types of edge iterators*/

        //! Beginning of the edges in the graph
        //! \return edge iterator pointing to the first edge in the graph
        edge_iterator edges_begin()
        {
            return edge_iterator(boost::edges(_graph));
        }

        //! Ending of the edges in the graph
        //! \return edge iterator pointing to the end of the edges in the graph
        edge_iterator edges_end()
        {
            return edge_iterator(boost::edges(_graph).second, boost::edges(_graph).second);
        }


        //! Beginning of the edges going out of a node in the graph
        //! \param[in] n : source node
        //! \return edge iterator pointing to the first out edge for the node
        out_edge_iterator out_edges_begin(Node n)
        {
            return out_edge_iterator(boost::out_edges(n, _graph));
        }

        //! Ending of the edges going out of a node in the graph
        //! \param[in] n : source node
        //! \return edge iterator pointing to the end of the out edges  for the node
        out_edge_iterator out_edges_end(Node n)
        {
            return out_edge_iterator(boost::out_edges(n, _graph).second, boost::out_edges(n, _graph).second);
        }

        //! Beginning of the edges going in of a node in the graph
        //! \param[in] n : target node
        //! \return edge iterator pointing to the first in edge for the node
        in_edge_iterator in_edges_begin(Node n)
        {
            return in_edge_iterator(boost::in_edges(n, _graph));
        }

        //! Ending of the edges going in of a node in the graph
        //! \param[in] n : target node
        //! \return edge iterator pointing to the end of the in edges  for the node
        in_edge_iterator in_edges_end(Node n)
        {
            return in_edge_iterator(boost::in_edges(n, _graph).second, boost::in_edges(n, _graph).second);
        }

    protected:
        //! Main graph structure
        GraphContainer _graph;
    };

} // ns graph
} // ns math
} // ns imageplus

#endif // HPP

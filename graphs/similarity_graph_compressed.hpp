// --------------------------------------------------------------
// Copyright (C)
// Universitat Politecnica de Catalunya (UPC) - Barcelona - Spain
// --------------------------------------------------------------

/*!
 * \file    similarity_graph.hpp
 * \author  Carles Ventura
 * \date    Apr 12, 2012
 *
 * \test    similarity_graph.test
 * \example similarity_graph.test
 *          A test-example for the SimilarityGraph class
 */

#ifndef IMAGEPLUS_MATH_GRAPH_SIMILARITY_GRAPH_HPP
#define IMAGEPLUS_MATH_GRAPH_SIMILARITY_GRAPH_HPP

#include <imageplus/math/graphs/graph.hpp>
#include <imageplus/math/graphs/algorithms/connected_components.hpp>

namespace imageplus
{
namespace math
{
namespace graphs
{
    //! Class for Node Properties
    template<class VDModel>
    class NodeDescProperties: public NodeDefaultProperties
    {
    public:
        //! Descriptor value
        VDModel value;
        //! Value for identifying the asset type (0 for intra asset ID and > 0 for inter assets IDs)
        uint64 asset_id;

        //! Default constructor
        NodeDescProperties()
        :   NodeDefaultProperties(),
            value(),
            asset_id(0)
        {}
    };

    //! Class for the Similarity Graph. Elements in the vertices are represented by a VDModel descriptor, which is
    //! used to compare elements according to the distance_functor.
    template<class VDModel, class distance_functor>
    class SimilarityGraph : public BoostGraph<kGraphUndirected, NodeDescProperties<VDModel>, EdgeDefaultProperties>
    {
    public:

        //! typedef for BoostGraph type
        typedef BoostGraph<kGraphUndirected, NodeDescProperties<VDModel>, EdgeDefaultProperties> GraphType;

        //! Method for inserting a new element in the similarity graph
        //!
        //! \param[in] desc: Descriptor value of the element to be inserted
        //! \param[in] name: String that identifies the element to be inserted
        //! \param[in] id: Label that identifies the element to be inserted
        //! \param[in] asset_id: ID for the asset which the element belongs to
        //!
        //! \returns the inserted Node
        typename GraphType::Node insert_node(VDModel& desc, std::string name, uint32 id, uint32 asset_id)
        {
            typename GraphType::Node node = this->add_node();
            this->node_properties(node).value = desc;
            this->node_properties(node).name = name;
            this->node_properties(node).id = id;
            this->node_properties(node).asset_id = asset_id;
            return node;
        }

        //! Method for computing the similarity graph. The similarity values are computed for each pair of elements.
        void calculate()
        {
            distance_functor dist;
            for(typename GraphType::node_iterator it_source = this->nodes_begin(); it_source != this->nodes_end(); ++it_source)
            {
                for(typename GraphType::node_iterator it_target = this->nodes_begin(); it_target != this->nodes_end(); ++it_target)
                {
                    if(this->node_properties(*it_source).id < this->node_properties(*it_target).id)
                    {
                        typename GraphType::EdgePropertiesType edge_props;
                        edge_props.weight = dist(this->node_properties(*it_source).value, this->node_properties(*it_target).value);
                        this->add_edge(*it_source, *it_target, edge_props);
                    }
                }
            }
        }

        //! Method for computing the similarity graph. The similarity values are computed for each pair of elements. Only edges between
        //! pairs of elements with similarity values greater than a threshold are considered.
        //!
        //! \param[in] threshold: If the similarity value between a pair of elements is below this threshold, then the edge that connects these two elements is not added.
        //!
        void calculate(float64 threshold)
        {
            distance_functor sim;
            for(typename GraphType::node_iterator it_source = this->nodes_begin(); it_source != this->nodes_end(); ++it_source)
            {
                for(typename GraphType::node_iterator it_target = this->nodes_begin(); it_target != this->nodes_end(); ++it_target)
                {
                    if(this->node_properties(*it_source).id < this->node_properties(*it_target).id)
                    {
                        float64 similarity = sim(this->node_properties(*it_source).value, this->node_properties(*it_target).value);
                        if (similarity > threshold)
                        {
                            typename GraphType::EdgePropertiesType edge_props;
                            edge_props.weight = similarity;
                            this->add_edge(*it_source, *it_target, edge_props);
                        }
                    }
                }
            }
        }


        //! Random walk algorithm that assigns each element a score depending on the relationship with
        //! the other elements in the graph. Random walk works as follows: Given a graph with vertices and a set of weighted
        //! edges, the ranking scores correspond to the likelihood of arriving in each of the vertices by traversing through the
        //! graph (with a random starting point), where the decision to take a particular path is defined by the weighted edges.
        //!
        //! \param[in] alpha: Weighting parameter in the random walk algorithm. The scores are biased by a vector of initial values with a factor (1-alpha).
        //!
        //! \return It returns a map with the score assigned to each node in the similarity graph
        std::map<typename GraphType::Node, float64> random_walk(float64 alpha = 0.8)
        {
            std::map<typename GraphType::Node, float64> normalization_values;
            for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
            {
                float64 sum_weights = 0;
                for(typename GraphType::node_iterator it2 = this->nodes_begin(); it2 != this->nodes_end(); ++it2)
                {
                    if(this->edge_exists(*it, *it2))
                    {
                        sum_weights += this->edge_properties(this->edge(*it,*it2)).weight;
                    }
                }
                normalization_values.insert(std::pair<typename GraphType::Node, float64>(*it, sum_weights));
            }

            uint64 size = this->num_nodes();
            std::map<typename GraphType::Node, float64> scores;
            std::map<typename GraphType::Node, float64> initial_scores;
            for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
            {
                scores.insert(std::pair<typename GraphType::Node, float64>(*it,1./size));
                initial_scores.insert(std::pair<typename GraphType::Node, float64>(*it,1./size));
            }
            uint64 num_iterations = 10;

            for(uint64 cur_iteration = 0; cur_iteration < num_iterations; cur_iteration++)
            {
                std::map<typename GraphType::Node, float64> new_scores;
                float64 sum_scores = 0;
                for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
                {
                    float64 cur_score = 0;
                    for(typename GraphType::node_iterator it2 = this->nodes_begin(); it2 != this->nodes_end(); ++it2)
                    {
                        if(this->edge_exists(*it, *it2))
                        {
                            float64 prob_trans = this->edge_properties(this->edge(*it,*it2)).weight;
                            if(normalization_values[*it2] != 0)
                            {
                                cur_score += prob_trans*scores[*it2]/normalization_values[*it2];
                            }
                            else //For directional graphs where a node can have no outgoing edges
                            {
                                cur_score += 1./size*scores[*it2];
                            }
                        }
                    }
                    cur_score = alpha*cur_score + (1-alpha)*initial_scores[*it];
                    sum_scores += cur_score;
                    new_scores.insert(std::pair<typename GraphType::Node, float64>(*it,cur_score));
                }

                //Normalizing new_scores
                for(typename std::map<typename GraphType::Node, float64>::iterator it = new_scores.begin(); it != new_scores.end(); ++it)
                {
                    new_scores[it->first] = it->second / sum_scores;
                }

                //Updating scores
                scores = new_scores;
            }
            return scores;
        }


        //! Random walk algorithm that assigns each element a score depending on the relationship with
        //! the other elements in the graph.
        //!
        //! \param[in] intra_weight: Weight for the elements belonging to the intra asset.
        //! \param[in] alpha: Weighting parameter in the random walk algorithm. The scores are biased by a vector of initial values with a factor (1-alpha).
        //!
        //! \return It returns a map with the score assigned to each node in the similarity graph
        std::map<typename GraphType::Node, float64> random_walk_keyframe_selection(float64 intra_weight = 0.5, float64 alpha = 0.8 )
        {
            uint64 intra_size = 0;
            std::map<typename GraphType::Node, float64> normalization_values;
            for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
            {
                float64 sum_weights = 0;
                if(this->node_properties(*it).asset_id == 1)
                {
                    intra_size++;
                }
                for(typename GraphType::node_iterator it2 = this->nodes_begin(); it2 != this->nodes_end(); ++it2)
                {
                    if(this->edge_exists(*it, *it2))
                    {
                        sum_weights += this->edge_properties(this->edge(*it,*it2)).weight;
                    }
                }
                normalization_values.insert(std::pair<typename GraphType::Node, float64>(*it, sum_weights));
            }

            uint64 size = this->num_nodes();
            std::map<typename GraphType::Node, float64> scores;
            std::map<typename GraphType::Node, float64> initial_scores;
            std::map<typename GraphType::Node,uint64 > rank_map;
            uint64 cur_node = 0;
            for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
            {
                if(this->node_properties(*it).asset_id == 1)
                {
                    scores.insert(std::pair<typename GraphType::Node, float64>(*it,intra_weight/intra_size));
                    initial_scores.insert(std::pair<typename GraphType::Node, float64>(*it,intra_weight/intra_size));
                    rank_map.insert(std::pair<typename GraphType::Node, uint64>(*it, cur_node));
                }
                else
                {
                    scores.insert(std::pair<typename GraphType::Node, float64>(*it,(1-intra_weight)/(size-intra_size)));
                    initial_scores.insert(std::pair<typename GraphType::Node, float64>(*it,(1-intra_weight)/(size-intra_size)));
                    rank_map.insert(std::pair<typename GraphType::Node, uint64>(*it, cur_node));
                }
                cur_node++;
            }
            uint64 max_iterations = 20;
            bool no_change = false;
            for(uint64 cur_iteration = 0; cur_iteration < max_iterations && !no_change; cur_iteration++)
            {
                std::map<typename GraphType::Node, float64> new_scores;
                std::map<float64, typename GraphType::Node> new_rank;
                float64 sum_scores = 0;
                for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
                {
                    float64 cur_score = 0;

                    for(typename GraphType::node_iterator it2 = this->nodes_begin(); it2 != this->nodes_end(); ++it2)
                    {
                        if(this->edge_exists(*it, *it2))
                        {
                            float64 prob_trans = this->edge_properties(this->edge(*it,*it2)).weight;
                            if(normalization_values[*it2] != 0)
                            {
                                cur_score += prob_trans*scores[*it2]/normalization_values[*it2];
                            }
                        }
                    }
                    cur_score = alpha*cur_score + (1-alpha)*initial_scores[*it];
                    sum_scores += cur_score;
                    new_scores.insert(std::pair<typename GraphType::Node, float64>(*it,cur_score));
                }

                //Normalizing new_scores
                for(typename std::map<typename GraphType::Node, float64>::iterator it = new_scores.begin(); it != new_scores.end(); ++it)
                {
                    new_scores[it->first] = it->second / sum_scores;
                    new_rank.insert(std::pair<float64, typename GraphType::Node>(it->second / sum_scores,it->first));
                }

                std::map<typename GraphType::Node,uint64 > new_rank_map;
                cur_node = 0;
                for(typename std::map<float64, typename GraphType::Node>::iterator it = new_rank.begin(); it != new_rank.end(); ++it)
                {
                    new_rank_map.insert(std::pair<typename GraphType::Node, uint64>(it->second, cur_node));
                    cur_node++;
                }


                no_change = true;
                for(typename std::map<typename GraphType::Node,uint64 >::iterator it = new_rank_map.begin(); it != new_rank_map.end(); ++it)
                {
                     if(it->second != rank_map[it->first])
                     {
                         no_change = false;
                     }
                }

                //Updating scores
                scores = new_scores;
                rank_map = new_rank_map;
            }

            return scores;
        }


        //! Mutual reinforcement algorithm also assigns each element a score depending on the relationship with
        //! the other elements in the graph. The greater the scores of the elements to which a node is connected and
        //! the greater the similarity values that connect them, the greater the score given to that node.
        //! There are two stages: (i) Mutual reinforcement-based score using only the subgraph containing the
        //! elements belonging to the intra video asset, and (ii) linear fusion taking into account the similarity
        //! to the elements belonging to the other video assets.
        //! More details can be found in: Ventura, C.; Giro-i-Nieto, X.; Vilaplana, V.; Giribet, D.; Carasusan, E., "Automatic keyframe selection based on mutual reinforcement algorithm."
        //!
        //! \param[in] intra_weight: Weight for the score achieved in the first stage(Mutual reinforcement-based intra video asset)
        //!
        //! \return It returns a map with the score assigned to each node in the similarity graph
        std::map<typename GraphType::Node, float64> reinforcement_based_ranking(float64 intra_weight = 0.5)
        {
            uint64 size = this->num_nodes();
            std::map<typename GraphType::Node, float64> scores;
            std::map<typename GraphType::Node, float64> initial_scores;
            std::map<typename GraphType::Node,uint64 > rank_map;
            uint64 cur_node = 0;
            uint64 intra_size=0;
            for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
            {
                if(this->node_properties(*it).asset_id == 1)
                {
                    intra_size++;
                }
            }
            for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
            {
                if(this->node_properties(*it).asset_id == 1)
                {
                    scores.insert(std::pair<typename GraphType::Node, float64>(*it,1./intra_size));
                    initial_scores.insert(std::pair<typename GraphType::Node, float64>(*it,1./intra_size));
                    rank_map.insert(std::pair<typename GraphType::Node, uint64>(*it, cur_node));
                    cur_node++;

                }
            }
            uint64 max_iterations = 20;
            bool no_change = false;
            for(uint64 cur_iteration = 0; cur_iteration < max_iterations && !no_change; cur_iteration++)
            {
                std::map<typename GraphType::Node, float64> new_scores;
                std::map<float64, typename GraphType::Node> new_rank;
                float64 sum_scores = 0;
                for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
                {
                    if(this->node_properties(*it).asset_id == 1)
                    {
                        float64 cur_score = 0;

                        for(typename GraphType::node_iterator it2 = this->nodes_begin(); it2 != this->nodes_end(); ++it2)
                        {
                            if(this->node_properties(*it2).asset_id == 1)
                            {
                                if(this->edge_exists(*it, *it2))
                                {
                                    float64 prob_trans = this->edge_properties(this->edge(*it,*it2)).weight;
                                    cur_score += prob_trans*scores[*it2];

                                }
                            }

                        }
                        sum_scores += cur_score;
                        new_scores.insert(std::pair<typename GraphType::Node, float64>(*it,cur_score));
                    }
                }

                if(sum_scores!=0)
                {
                    //Normalizing new_scores
                    for(typename std::map<typename GraphType::Node, float64>::iterator it = new_scores.begin(); it != new_scores.end(); ++it)
                    {
                        new_scores[it->first] = it->second / sum_scores;
                        new_rank.insert(std::pair<float64, typename GraphType::Node>(it->second / sum_scores,it->first));
                    }

                    std::map<typename GraphType::Node,uint64 > new_rank_map;
                    cur_node = 0;
                    for(typename std::map<float64, typename GraphType::Node>::iterator it = new_rank.begin(); it != new_rank.end(); ++it)
                    {
                        new_rank_map.insert(std::pair<typename GraphType::Node, uint64>(it->second, cur_node));
                        cur_node++;
                    }


                    no_change = true;
                    for(typename std::map<typename GraphType::Node,uint64 >::iterator it = new_rank_map.begin(); it != new_rank_map.end(); ++it)
                    {
                         if(it->second != rank_map[it->first])
                         {
                             no_change = false;
                         }
                    }

                    //Updating scores
                    scores = new_scores;
                    rank_map = new_rank_map;
                }
                else
                {
                    scores = initial_scores;
                    no_change=true;
                }
            }


            //Linear fusion with external keyframes
            uint64 num_external_keyframes = size-intra_size;
            if(num_external_keyframes > 0)
            {
                float64 weighting_parameter = (1-intra_weight)/num_external_keyframes;
                float64 sum_scores = 0;
                for(typename GraphType::node_iterator it = this->nodes_begin(); it != this->nodes_end(); ++it)
                {
                    if(this->node_properties(*it).asset_id == 1)
                    {
                        float64 cur_score = intra_weight*scores[*it];
                        for(typename GraphType::node_iterator it2 = this->nodes_begin(); it2 != this->nodes_end(); ++it2)
                        {
                            if(this->node_properties(*it2).asset_id > 1)
                            {
                                float64 prob_trans = this->edge_properties(this->edge(*it,*it2)).weight;
                                cur_score += weighting_parameter*prob_trans;
                            }
                        }
                        scores[*it] = cur_score;
                        sum_scores += cur_score;
                    }
                }
                //Normalizing scores
                for(typename std::map<typename GraphType::Node, float64>::iterator it = scores.begin(); it != scores.end(); ++it)
                {
                    scores[it->first] = it->second / sum_scores;
                }
            }


            return scores;

        }

    };

} // ns graph
} // ns math
} // ns imageplus


#endif // HPP

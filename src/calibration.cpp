#include "msgpack.hpp"
#include "Rcpp.h"
#include "vector"
#include "set"
#include "map"
#include "string"
#include "iostream"
#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppMsgPack)]]
// [[Rcpp::depends(BH)]]

struct Node {
  Node * parent;
  int level;
  // string name;
  //std::vector<node_feature> features;
  std::set<int> active_idx;
  std::set<int> terminal_idx;
  std::map<std::string, Node*> children;
};


//TO DO: write out function signatures
// struct Node* NewNode(Node *parent, string name, set<int> active_idx); //Function for node creation
// struct Node* updateKeyWordTree(vector<string> clone_sequences, vector<int> ndn_start_vec, vector<int> ndn_end_vec, vector<unsigned char> existing_tree);
//Function to update the tree with a new sequence
//packTree // Function for serializing tree, so that it can be ported to R environment
//unpackTree //Function for de-serializing tree from R
//findMatchingIdx //Function to find compatible indices based on existing keyword tree

struct Node* NewNode(Node *parent, std::string name, std::set<int> active_idx, std::set<int> terminal_idx) {
  Node *node = new Node;
  node->active_idx = active_idx;
  node->terminal_idx = terminal_idx;
  if(parent == NULL) {
    node->level = 0;
  } else {
    node->level = parent->level + 1;
    parent->children.insert(std::pair<std::string, Node*>(name, node));
  }
  return node;
}

// packTree and unpackTree are serializations method so the tree can be saved into R or to disk if desired
void packTree(Node * node,  msgpack::packer<std::stringstream>& pkr) {
  pkr.pack_array(3); //active_idx, children

  std::set<int> active_idx = node->active_idx;
  std::vector<int> va(active_idx.begin(), active_idx.end());
  pkr.pack(va);
  
  std::set<int> terminal_idx = node->terminal_idx;
  std::vector<int> vt(terminal_idx.begin(), terminal_idx.end());
  pkr.pack(vt);
  
  std::map<std::string, Node*> children = node->children;
  pkr.pack_map(children.size());
  for (auto const& x : children) {
    pkr.pack(x.first);  // string key
    packTree(x.second, pkr); //node value
  }
}

Node* unpackTree(msgpack::object & node_obj, std::string name = "root", Node* parent = NULL) {
  std::vector<msgpack::object> node_obj_array;
  node_obj.convert(node_obj_array);
  
  //active_idx
  std::vector<int> v;
  node_obj_array[0].convert(v);
  std::set<int> active_idx(v.begin(), v.end());
  
  node_obj_array[1].convert(v);
  std::set<int> terminal_idx(v.begin(), v.end());
  
  Node* current_node = NewNode(parent, name, active_idx, terminal_idx);
  
  //children
  std::map<std::string, Node*> children;
  msgpack::object_kv* p = node_obj_array[2].via.map.ptr;
  msgpack::object_kv* const pend = node_obj_array[2].via.map.ptr + node_obj_array[2].via.map.size;
  for (; p < pend; ++p) {
    std::string name = p->key.as<std::string>();
    msgpack::object node_child_obj = p->val.as<msgpack::object>();
    Node* node_child = unpackTree(node_child_obj, name, current_node);
    children[name] = node_child;
  }
  return(current_node);
}


// [[Rcpp::export]]
Rcpp::RawVector C_updatekeyWordTree(std::vector<std::string> clone_sequences, std::vector<int> ndn_start_vec, std::vector<int> ndn_end_vec, std::vector<int> clone_idx, std::vector<unsigned char> packed_tree) {
  Node* root;
  if(packed_tree.size() == 0) { // E.g. raw vector of length zero from R
    root = NewNode(NULL, "root", {}, {});
  } else {
    std::string message(packed_tree.begin(), packed_tree.end());
    msgpack::object_handle oh;
    msgpack::unpack(oh, message.data(), message.size());
    msgpack::object obj = oh.get();
    root = unpackTree(obj, "root", NULL);
  }
  
  for(uint i=0; i<clone_sequences.size(); i++) {
    Node * current_node = root;
    for(uint j=0; j < clone_sequences[i].size(); j++) {
      char nucleotide = clone_sequences[i][j];
      std::string is_ndn = ((j >= ndn_start_vec[i] && j <= ndn_end_vec[i])) ? "1" : "0";
      std::string name = nucleotide + is_ndn;

      //test whether nucleotide+is_ndn child node exists
      if(current_node->children.count(name) == 1) {
        current_node = current_node->children[name]; // change current node to matching child node
        if(j == clone_sequences[i].size() - 1) {
          current_node->terminal_idx.insert(clone_idx[i]); // add terminal_idx at last nucleotide
        } else {
          current_node->active_idx.insert(clone_idx[i]); // add active_idx to matching child node
        }
      } else {
        if(j == clone_sequences[i].size() - 1) {
          current_node = NewNode(current_node, name, {}, {clone_idx[i]}); // create a new child node with terminal_idx and assign it to current node
        } else {
          current_node = NewNode(current_node, name, {clone_idx[i]}, {}); // create a new child node with active_idx and assign it to current node
        }
      }
    }
  }
  
  std::stringstream buffer;
  msgpack::packer<std::stringstream> pk(&buffer);
  packTree(root, pk);
  std::string bufstr = buffer.str();
  Rcpp::RawVector rawbuffer(bufstr.begin(), bufstr.end());
  return rawbuffer;
}

//Pre-order traversal with pruning when cost > max_dist
//result_idx is the output
void findSeq(std::string clone_sequence, Node* node, std::string name, double parent_cost, double max_dist, double ndn_cost, std::map<int, double> & result_idx) {
  double mismatch_cost = (name[1] == '1') ? ndn_cost : 1;
  double node_cost;
  if(name == "root") {
    node_cost = parent_cost;
  } else {
    node_cost = (name[0] == clone_sequence[node->level - 1]) ? parent_cost : (parent_cost + mismatch_cost);
  }

  if(node_cost <= max_dist) {
    //case 1 - if end of query
    if(node->level == clone_sequence.size()) {
      for (const int &idx : node->active_idx)
        result_idx[idx] = node_cost;
    }
    //case 2 - Any terminal_idx, e.g., is leaf
    if(node->terminal_idx.size() != 0) { // if(node->children.size() == 0)
      for (const int &idx : node->terminal_idx)
        result_idx[idx] = node_cost;
    }
    //case X - if query is not finished, recursively search children
    if(node->level <= clone_sequence.size()) {
      std::map<std::string, Node*> children = node->children;
      for (auto const& x : children)
        findSeq(clone_sequence, x.second, x.first, node_cost, max_dist, ndn_cost, result_idx);
    }
  }
}

// [[Rcpp::export]]
Rcpp::List C_findMatchingIdx(std::vector<std::string> clone_sequences, std::vector<double> max_dist, double ndn_cost, std::vector<unsigned char> packed_tree) {
  std::string message(packed_tree.begin(), packed_tree.end());
  msgpack::object_handle oh;
  msgpack::unpack(oh, message.data(), message.size());
  msgpack::object obj = oh.get();
  Node* root = unpackTree(obj, "root", NULL);
  
  Rcpp::List all_results(clone_sequences.size());
  for(uint i=0; i<clone_sequences.size(); i++) {
    std::map<int, double> result_idx;
    findSeq(clone_sequences[i], root, "root", 0, max_dist[i], ndn_cost, result_idx);
    // Rcpp::IntegerVector ri(result_idx.begin(), result_idx.end());
    // all_results[i] = ri;
    Rcpp::NumericVector all_results_cost(result_idx.size());
    Rcpp::CharacterVector all_results_idx(result_idx.size());
    int j = 0;
    for (auto const& x : result_idx) {
      all_results_idx[j] = std::to_string(x.first);
      all_results_cost[j] = x.second;
      j++;
    }
    all_results_cost.names() = all_results_idx;
    all_results[i] = all_results_cost;
  }
  
  return(all_results);
}

// Resolve familyTies
std::string resolveFamilyLocus(std::string& familyTies) {
  std::vector<std::string> split;
  boost::split(split, familyTies, boost::is_any_of(","));
  // check that each split is at least 3 characters, so we don't run into a substr error
  for(uint j=0;j<split.size();j++) {
    if(split[j].size() < 3) return("");
  }
  std::string locus_temp = split[0].substr(0,3);
  for(uint j=1;j<split.size();j++) {
    if(split[j].substr(0,3) != locus_temp) return("");
  }
  return(locus_temp);
}

// utility function for finding consensus locus
// [[Rcpp::export]]
std::vector<std::string> C_getLocus(std::vector<std::string> cloneResolved, std::vector<std::string> vFamilyName,
                               std::vector<std::string> vFamilyTies, std::vector<std::string> dFamilyName,
                               std::vector<std::string> dFamilyTies, std::vector<std::string> jFamilyName,
                               std::vector<std::string> jFamilyTies) {
  std::vector<std::string> locus(cloneResolved.size());
  for(uint i=0; i<cloneResolved.size(); i++) {
    locus[i] = ""; // default value
    // First check for annotations in primary labels
    if(vFamilyName[i].size() > 0 & vFamilyName[i].find("IG") == 0) {
      locus[i] = vFamilyName[i].substr(0,3);
      continue;
    } else if(dFamilyName[i].size() > 0 & dFamilyName[i].find("IG") == 0) {
      locus[i] = dFamilyName[i].substr(0,3);
      continue;
    } else if(jFamilyName[i].size() > 0 & jFamilyName[i].find("IG") == 0) {
      locus[i] = jFamilyName[i].substr(0,3);
      continue;
    } else {
      if(vFamilyTies[i].size() > 0) {
        locus[i] = resolveFamilyLocus(vFamilyTies[i]);
        if(locus[i] != "") continue;
      }
      if(dFamilyTies[i].size() > 0) {
        locus[i] = resolveFamilyLocus(dFamilyTies[i]);
        if(locus[i] != "") continue;
      }
      if(jFamilyTies[i].size() > 0) {
        locus[i] = resolveFamilyLocus(jFamilyTies[i]);
      }
    }
  }
  // What is the (biological) difference between IGH and IGH_D?
  for(uint i=0; i<cloneResolved.size(); i++) {
    if(locus[i] == "IGH" & cloneResolved[i] == "DJ") locus[i] = "IGH_D";
  }
  return(locus);
}

// utility function for calculating cdr3
// [[Rcpp::export]]
std::vector<int> C_determineCdr3Index(std::vector<int> cloneLength, std::vector<int> vIndex, std::vector<int> dIndex, std::vector<int> cdr3Length) {
  std::vector<int> cdr3Index(cloneLength.size());
  for(uint i=0; i<cloneLength.size(); i++) {
    cdr3Index[i] = vIndex[i];
    if(cdr3Index[i] == -1) cdr3Index[i] = dIndex[i];
    if(cdr3Index[i] != -1) cdr3Index[i] += cdr3Length[i];
    if(cdr3Index[i] < 0 | cdr3Index[i] >= cloneLength[i]) cdr3Index[i] = -1;
  }
  return(cdr3Index);
}


// void test() {
//   Node * root = NewNode(NULL, "root", {});
//   Node * child = NewNode(root, "B1", {1,2,3});
//   cerr << root->level << "\n";
//   cerr << child->level << "\n";
//   cerr << root->children["B1"]->level << "\n";
// }

// void print(vector<unsigned char> packed_tree) {
//   Node * root = NewNode(NULL, "root", {});
//   Node * child = NewNode(root, "B1", {1,2,3});
//   cerr << root->level << "\n";
//   cerr << child->level << "\n";
//   cerr << root->children["B1"]->level << "\n";
// }




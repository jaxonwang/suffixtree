#include <array>
#include <charconv>
#include <sstream>
#include <chrono>
#include <iostream>
#include <memory>
#include <unordered_map>
using namespace std;
using namespace std::chrono;

class SuffixTree;
class Node;
class NonLeaf;
class InternalNode;
class Leaf;
class Edge;
class InternalEdge;
class LeafEdge;

using CharT = char;

class Node {
public:
  using edge_ptr = std::unique_ptr<Edge>;
  using node_ptr = std::shared_ptr<Node>;
  Node() = default;
  virtual ~Node() = default;
  // virtual node_ptr next_node(const CharT c) = 0;
  // virtual edge_ptr &next_edge(const CharT c) = 0;
};

class NonLeaf : public Node {
public:
  using edge_ptr = std::unique_ptr<Edge>;
  using node_ptr = std::shared_ptr<Node>;
  std::unordered_map<CharT, node_ptr> children;
  std::unordered_map<CharT, edge_ptr> adj_edges;
  node_ptr next_node(const CharT c) {
    if (children.count(c))
      return children[c];
    return nullptr;
  }
  edge_ptr &next_edge(const CharT c) {
    static edge_ptr empty_edge{};
    if (adj_edges.count(c))
      return adj_edges[c];
    return empty_edge;
  }
};

class InternalNode : public NonLeaf {
  friend class SuffixTree;
  using internal_node_ptr = std::shared_ptr<InternalNode>;
  using node_ptr = std::shared_ptr<Node>;
  using suffix_link_t = std::weak_ptr<InternalNode>;
  suffix_link_t suffix_link;

public:
  InternalNode() = default;
  void set_suffix_link(const internal_node_ptr &p) {
    suffix_link = suffix_link_t(p);
  }
  suffix_link_t get_suffix_link() { return suffix_link; }
};

class Leaf : public Node {
  using size_type = size_t;

public:
  size_type length;
  Leaf() : length(0) {}
};

class Edge {
public:
  using edge_ptr = std::unique_ptr<Edge>;
  using size_type = size_t;

  Edge() = default;
  virtual ~Edge() = default;
  virtual size_type start_position() = 0;
  virtual size_type end_position() = 0;
  virtual size_type length() { return end_position() - start_position() + 1; }
  virtual std::pair<edge_ptr, edge_ptr>
  break_edge(const size_type pos) const = 0;
};

class InternalEdge : public Edge {
  using size_type = typename Edge::size_type;
  using edge_ptr = typename Edge::edge_ptr;
  const size_type start_pos;
  const size_type end_pos;

public:
  InternalEdge(const size_type start_pos, const size_type end_pos)
      : start_pos(start_pos), end_pos(end_pos) {}
  size_type start_position() override { return start_pos; }
  size_type end_position() override { return end_pos; }
  virtual std::pair<edge_ptr, edge_ptr>
  break_edge(const size_type pos) const override {
    if (pos < start_pos || pos >= end_pos)
      throw std::logic_error("bad pos");
    auto new_first = edge_ptr{new InternalEdge(start_pos, pos)};
    auto new_second = edge_ptr{new InternalEdge(pos + 1, end_pos)};
    return make_pair(move(new_first), move(new_second));
  }
};

class LeafEdge : public Edge {
  using size_type = typename Edge::size_type;
  using edge_ptr = typename Edge::edge_ptr;
  const size_type start_pos;
  const size_type &end_pos;

public:
  LeafEdge(const size_type start_pos, const size_type &end_pos)
      : start_pos(start_pos), end_pos(end_pos) {}
  virtual size_type start_position() override { return start_pos; }
  virtual size_type end_position() override { return end_pos; }
  virtual std::pair<edge_ptr, edge_ptr>
  break_edge(const size_type pos) const override {
    if (pos < start_pos || pos >= end_pos)
      throw std::logic_error("bad pos");
    auto new_first = edge_ptr{new InternalEdge(start_pos, pos)};
    auto new_second = edge_ptr{new LeafEdge(pos + 1, end_pos)};
    return make_pair(move(new_first), move(new_second));
  }
};

class SuffixTreeBase {
public:
  using edge_ptr = std::unique_ptr<Edge>;
  using node_ptr = std::shared_ptr<Node>;
  using internal_node_ptr = std::shared_ptr<InternalNode>;
  using size_type = size_t;

protected:
  using StringT = std::string;
  const StringT text;
  const size_type textlen;
  static const CharT term_symbol = '\0';
  internal_node_ptr root;

public:
  SuffixTreeBase() = delete;
  SuffixTreeBase(const StringT &text)
      : text(text), textlen(text.size()),
        root(std::make_shared<InternalNode>()) {
    root->set_suffix_link(root);
  }
  internal_node_ptr break_edge(internal_node_ptr &edge_owner, const CharT c,
                               const size_type pos) {
    // break edge[start, end] into [start, pos] [post+1, end]
    auto real_node = edge_owner;
    auto &e = real_node->adj_edges[c];
    auto edge_pair = e->break_edge(pos);
    auto first_edge = move(edge_pair.first);
    auto second_edge = move(edge_pair.second);
    auto new_internal = std::make_shared<InternalNode>();

    const CharT breakpoint_c = text[pos + 1];
    new_internal->children[breakpoint_c] = real_node->children[c];
    new_internal->adj_edges[breakpoint_c] = move(second_edge);

    real_node->children[c] = std::static_pointer_cast<Node>(
        std::shared_ptr<InternalNode>{new_internal});
    real_node->adj_edges[c] = move(first_edge);
    return new_internal;
  }
  void new_leaf(internal_node_ptr &current_node, const size_type start_pos,
                const size_type &end_pos) {
    const CharT c = text[start_pos];
    current_node->adj_edges[c] =
        std::unique_ptr<Edge>{new LeafEdge{start_pos, end_pos}};
    current_node->children[c] =
        std::static_pointer_cast<Node>(std::make_shared<Leaf>());
  }

  void _print(const internal_node_ptr &node) const {
    for (auto &kv : node->children) {
      CharT c = kv.first;
      std::cout << node.get() << " ";
      auto &edge = node->adj_edges[c];
      size_type end_pos = edge->end_position();
      for (size_type i = edge->start_position(); i <= end_pos; i++) {
        if (text[i] == term_symbol)
          std::cout << '$';
        else
          std::cout << text[i];
      }
      std::cout << " " << kv.second.get() << std::endl;
    }
    internal_node_ptr suffix_node = node->get_suffix_link().lock();
    if (suffix_node) {
      std::cout << node.get() << " -> " << suffix_node.get() << endl;
    }
    for (auto &kv : node->children) {
      auto next = std::dynamic_pointer_cast<InternalNode>(kv.second);
      if (next)
        _print(next);
    }
  }

  bool static _same_tree(const internal_node_ptr &root1,
                         const internal_node_ptr &root2) {
    if (root1->adj_edges.size() != root2->adj_edges.size()) {
      std::cout << "size different" << std::endl;
      return false;
    }
    for (auto &c : root1->adj_edges) {
      if (root2->adj_edges.count(c.first) == 0) {
        std::cout << c.first << " not in the second" << std::endl;
        return false;
      }
      auto &edge1 = c.second;
      auto &edge2 = root2->adj_edges[c.first];
      if (edge1->start_position() != edge2->start_position()) {
        std::cout << "Different start position" << std::endl;
        return false;
      }
      if (edge1->end_position() != edge2->end_position()) {
        std::cout << "Different end position" << std::endl;
        return false;
      }
    }
    for (auto &c : root1->children) {
      if (root2->children.count(c.first) == 0) {
        std::cout << "children size different" << std::endl;
        return false;
      }
      auto child1 = std::dynamic_pointer_cast<InternalNode>(c.second);
      auto child2 =
          std::dynamic_pointer_cast<InternalNode>(root2->children[c.first]);
      // if child is leaf then the cooresponding child should be leaf
      if ((!child1 && child2) || (child1 && !child2)) {
        std::cout << "different child type" << std::endl;
        return false;
      }
      if (!child1)
        continue;
      if (!_same_tree(child1, child2))
        return false;
    }
    return true;
  }

  bool operator==(const SuffixTreeBase &other) const {
    return _same_tree(root, other.root);
  }

  void print() const { _print(root); }
};

class SuffixTreeNaive : public SuffixTreeBase {

public:
  SuffixTreeNaive(const StringT &s) : SuffixTreeBase(s) {}
  void build() {
    for (size_t i = 0; i <= textlen; i++) {
      int j = i;
      internal_node_ptr current_node = root;

      while (true) { // will eventually create a leaf
        CharT current_char = text[j];

        auto &nextedge_ptr = current_node->next_edge(current_char);
        if (!nextedge_ptr) {
          new_leaf(current_node, j, textlen);
          goto end_insert_i;
        }

        size_t start = nextedge_ptr->start_position();
        size_t end = nextedge_ptr->end_position();
        for (; start <= end; start++, j++) {
          if (text[start] != text[j]) {
            auto new_node = break_edge(current_node, current_char, start - 1);
            new_leaf(new_node, j, textlen); // include the term symbol
            goto end_insert_i;
          }
        }
        current_node = std::dynamic_pointer_cast<InternalNode>(
            current_node->next_node(current_char));
        if (!current_node)
          throw std::logic_error("Goes to Leaf!");
      }

    end_insert_i:;
    }
  }
};

class SuffixTreeImpl : public SuffixTreeBase {
private:
  enum class ExtensionType : int {
    newleaf,
    newnode,
    extendleaf, // extend the end point of leaf or just do nothing for implicit
                // suffix tree
  };
  size_type leaf_end_pos;

public:
  SuffixTreeImpl(const StringT &s) : SuffixTreeBase(s), leaf_end_pos(0) {}

  // try travel down to the point S[j, i)(could be empty) is matchd, and append
  // S(i), create node if necessary
  ExtensionType travel_down_and_extend(internal_node_ptr &node,
                                       const size_type _j, const size_type i,
                                       internal_node_ptr &stop_node,
                                       internal_node_ptr &suffix_link_node,
                                       size_type &travel_depth) {
    size_type j = _j;
    auto current_node = node;
    ExtensionType ret = ExtensionType::extendleaf;
    suffix_link_node = nullptr;
    travel_depth = 0;

    for (;;) {
      const CharT c = text[j];
      auto &edge = current_node->next_edge(c);
      if (!edge) { // node has no edge with label start with c
        new_leaf(current_node, j, leaf_end_pos);
        ret = ExtensionType::newleaf;
        break;
      }

      size_type remain_length = i - j + 1;
      size_type start_pos = edge->start_position();
      size_type edge_length = edge->length();
      if (remain_length <= edge_length) { // travel will end in this edge
        size_type possible_diverge_point = start_pos + remain_length - 1;
        if (text[possible_diverge_point] == text[i]) {
          break; // nothing to do, break and will return extendleaf
        } else {
          auto new_node =
              break_edge(current_node, c, possible_diverge_point - 1);
          new_leaf(new_node, i, leaf_end_pos);
          ret = ExtensionType::newnode;
          suffix_link_node = new_node;
          break;
        }
      }
      // skip this node and travel to next node
      j += edge_length;
      travel_depth += edge_length;
      auto maybe_internal =
          std::dynamic_pointer_cast<InternalNode>(current_node->next_node(c));
      if (!maybe_internal) { // extend on leaf edge
        if (i - j + 1 != 1)
          throw("remain error");
        break;
      } else {
        current_node = maybe_internal; // go to next node
      }
    }
    stop_node = current_node;
    if (!suffix_link_node)
      suffix_link_node = stop_node;
    return ret;
  }

  void build() { // S[j,i)
    internal_node_ptr longest_node = root;
    internal_node_ptr node_to_start = longest_node;

    size_type j_start_with = 0;
    size_type start_depth = 0;

    for (size_type i = 0; i <= textlen; i++) {
      leaf_end_pos = i;
      internal_node_ptr node_to_update_suffix_link = nullptr;

      // find the depth of the longest_node

      for (size_type j = j_start_with; j <= i; j++) {

        internal_node_ptr stop_node;
        internal_node_ptr suffix_link_node;
        size_type travel_depth;
        auto ret =
            travel_down_and_extend(node_to_start, j + start_depth, i, stop_node,
                                   suffix_link_node, travel_depth);

        if (node_to_update_suffix_link) {
          // set suffix link
          node_to_update_suffix_link->set_suffix_link(suffix_link_node);
        }

        if (ret == ExtensionType::newnode) { // need to specify suffix link must
                                             // be updated
          node_to_update_suffix_link = suffix_link_node;
        } else {
          node_to_update_suffix_link = nullptr;
        }

        if (ret != ExtensionType::extendleaf) { // new leaf is created
          j_start_with = j + 1;
          node_to_start = stop_node->get_suffix_link().lock(); // next to start
          start_depth += (stop_node.get() == root.get()) ? 0 : travel_depth - 1;
          cout << start_depth << endl;
        } else { // no need for the future extend , since they must be extesion
                 // type 3
          break;
        }
      }
    }
  }
};

string rand_str(size_t length){
   stringstream ss; 
   for (size_t i = 0; i < length; i++) {
    ss << static_cast<char>('a' + rand() % 26);
   }
   return ss.str();
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cout << "Usage: a.out str" << std::endl;
    return -1;
  }

  std::string s1 = rand_str(stoi(argv[1]));

  auto t_start = steady_clock::now();
  SuffixTreeNaive st2{s1};
  st2.build();

  auto t_end = steady_clock::now();
  auto dur = duration<double>(t_end - t_start).count();
  cout << "Test: naive Execution time: " << dur << endl;

  t_start = steady_clock::now();
  SuffixTreeImpl st1{s1};
  st1.build();
  t_end = steady_clock::now();
  dur = duration<double>(t_end - t_start).count();
  cout << "Test: impl Execution time: " << dur << endl;

  std::cout << (st2 == st1) << std::endl;

  return 0;
}

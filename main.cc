#include <iostream>
#include <memory>
#include <unordered_map>

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
  using node_ptr = std::shared_ptr<Node>;
  using suffix_link_t = std::weak_ptr<Node>;
  suffix_link_t suffix_link;

public:
  InternalNode() = default;
  void set_suffix_link(const node_ptr &p) { suffix_link = suffix_link_t(p); }
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
  virtual size_type length() { return end_position() - start_position(); }
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
  static const CharT term_symbol = '\0';
  internal_node_ptr root;

public:
  SuffixTreeBase() = delete;
  SuffixTreeBase(const StringT &text) : text(text) {}
  internal_node_ptr break_edge(internal_node_ptr &edge_owner, const CharT c,
                               const size_type pos) {
    auto real_node = edge_owner;
    auto &e = real_node->adj_edges[c];
    auto edge_pair = e->break_edge(pos);
    auto first_edge = move(edge_pair.first);
    auto second_edge = move(edge_pair.second);
    auto new_internal = std::unique_ptr<InternalNode>{new InternalNode{}};

    new_internal->children[c] = real_node->children[c];
    new_internal->adj_edges[c] = move(second_edge);

    real_node->children[c] = std::static_pointer_cast<Node>(
        std::shared_ptr<InternalNode>{move(new_internal)});
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
    for (auto &kv : node->children) {
      auto next = std::dynamic_pointer_cast<InternalNode>(kv.second);
      if (next)
        _print(next);
    }
  }

  void print() const { _print(root); }
};

class SuffixTreeNaive : public SuffixTreeBase {
  size_t textlen;

public:
  SuffixTreeNaive(const StringT &s) : SuffixTreeBase(s) {}
  void build() {
    textlen = text.size();
    root = std::make_shared<InternalNode>();
    for (size_t i = 0; i < textlen; i++) {
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
        size_t end = nextedge_ptr->start_position();
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

int main() {

  std::string s1{"banana"};
  SuffixTreeNaive st1{s1};
  st1.build();
  st1.print();

  return 0;
}

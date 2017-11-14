#define BOOST_PYTHON_STATIC_LIB

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <utility>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "matrix/matrix.hpp"

namespace py = boost::python;

namespace seqbdd {
    // types
    typedef struct _Node {
        char label;
        const struct _Node* branch1;
        const struct _Node* branch0;

        _Node(char l, const struct _Node* b1, const struct _Node* b0)
            : label(l), branch1(b1), branch0(b0) {}
    } Node;

    typedef std::tuple<char, const Node*, const Node*> NodeKey;

    class nodehash {
    public:
        size_t operator()(const NodeKey& x) const {
            return (size_t)std::get<1>(x) ^ (size_t)std::get<2>(x);
        }
    };

    // static variables
    static const Node TERM0('0', NULL, NULL);
    static const Node TERM1('1', NULL, NULL);
    static std::unordered_map<NodeKey, std::unique_ptr<const Node>, nodehash> cache;

    // functions
    const Node* get_node(char x, const Node* p1, const Node* p0) {
        if (p1 == &TERM0) {
            return p0;
        } else {
            NodeKey key(x, p1, p0);

            if (!cache.count(key)) {
                cache[key] = std::unique_ptr<const Node>(new Node(x, p1, p0));
            }

            return cache[key].get();
        }
    }

    void _get_sequences_trace(const Node* node, std::vector<std::string>& sequences,
                          std::string& queue) {
        if (node == &TERM1) {
            sequences.push_back(queue);
        } else if (node != &TERM0) {
            _get_sequences_trace(node->branch0, sequences, queue);

            queue.push_back(node->label);
            _get_sequences_trace(node->branch1, sequences, queue);
            queue.pop_back();
        }
    }

    std::vector<std::string> get_sequences(const Node* root) {
        std::vector<std::string> sequences;
        std::string queue;

        _get_sequences_trace(root, sequences, queue);

        return sequences;
    }

    inline int sign(const Node* x) {
        return (x == &TERM0)? 0 : 1;
    }

    template<typename Func>
    const Node* meld(Func op, const Node* p, const Node* q) {
        if (p == &TERM0 || q == &TERM0 || p == q) {
            if (!op(sign(p), sign(q))) {
                return &TERM0;
            } else if (p != &TERM0) {
                return p;
            } else {
                return q;
            }
        }

        char x = p->label;
        char y = q->label;

        if (x > y) {
            return get_node(x, meld(op, p->branch1, &TERM0), meld(op, p->branch0, q));
        } else if (x < y) {
            return get_node(y, meld(op, &TERM0, q->branch1), meld(op, p, q->branch0));
        } else {
            return get_node(x, meld(op, p->branch1, q->branch1), meld(op, p->branch0, q->branch0));
        }
    }

    const Node* union_(const Node* p, const Node* q) {
        return meld([](int x, int y) -> int { return x+y; }, p, q);
    }

    const Node* intersection(const Node* p, const Node* q) {
        return meld([](int x, int y) -> int { return x*y; }, p, q);
    }

    const Node* difference(const Node* p, const Node* q) {
        return meld([](int x, int y) -> int { return x&(!y); }, p, q);
    }

    const Node* symmetric_difference(const Node* p, const Node* q) {
        return meld([](int x, int y) -> int { return x^y; }, p, q);
    }

    const Node* one_terminal() {
        return &TERM1;
    }

    const Node* zero_terminal() {
        return &TERM0;
    }

    const Node* onset(const Node* root, char x) {
        const Node* node = root;
        while (node != &TERM1 && node != &TERM0) {
            if (node->label == x) {
                return node->branch1;
            }
            node = node->branch0;
        }

        return node;
    }

    const Node* push(const Node* root, char x) {
        return get_node(x, root, &TERM0);
    }

    char top(const Node* root) {
        return root->label;
    }

    int count(const Node* root) {
        if (root == &TERM0) {
            return 0;
        } else if (root == &TERM1) {
            return 1;
        }

        const Node* node = root;
        std::unordered_map<const Node*, int> counter { { &TERM0, 0 }, { &TERM1, 1 } };
        std::vector<const Node*> stack { root };
        while (!stack.empty()) {
            node = stack.back();

            if (counter.count(node->branch0) == 0) {
                stack.push_back(node->branch0);
            } else if (counter.count(node->branch1) == 0) {
                stack.push_back(node->branch1);
            } else {
                counter[node] = counter[node->branch0] + counter[node->branch1];
                stack.pop_back();
            }
        }

        return counter[root];
    }

    int count_node(const Node* root) {
        if (root == &TERM0 || root == &TERM1) {
            return 0;
        }

        const Node* node = root;
        std::unordered_set<const Node*> set { &TERM0, &TERM1 };
        std::vector<const Node*> stack { root };
        while (!stack.empty()) {
            node = stack.back();

            if (set.count(node->branch0) == 0) {
                stack.push_back(node->branch0);
            } else if (set.count(node->branch1) == 0) {
                stack.push_back(node->branch1);
            } else {
                set.insert(node);
                stack.pop_back();
            }
        }

        return set.size() - 2;
    }

    bool equal_nodes(const Node* root, const Node* other) {
        return root == other;
    }

    bool not_equal_nodes(const Node* root, const Node* other) {
        return root != other;
    }

    const Node* make_string(const py::str& sequence) {
        const Node* root = &TERM1;
        for (int i=py::len(sequence)-1; i >= 0; --i) {
            root = get_node(py::extract<char>(sequence[i]), root, &TERM0);
        }

        return root;
    }

    const Node* from_sequences(const py::list& sequences) {
        const Node* root = &TERM0;
        for (int i=0; i < py::len(sequences); ++i) {
            root = union_(root, make_string(py::extract<py::str>(sequences[i])));
        }

        return root;
    }

    void _values_trace(const Node* node, py::list& sequences, std::string& queue) {
        if (node == &TERM1) {
            sequences.append(py::str(std::string(queue.begin(), queue.end()).c_str()));
        } else if (node != &TERM0) {
            _values_trace(node->branch0, sequences, queue);

            queue.push_back(node->label);
            _values_trace(node->branch1, sequences, queue);
            queue.pop_back();
        }
    }

    py::list values(const Node* root) {
        py::list sequences;
        std::string queue;

        _values_trace(root, sequences, queue);

        return sequences;
    }

    bool has_sequence(const Node* root, const py::str& sequence) {
        bool match = false;

        const Node* node = root;
        for (int i=0; i < py::len(sequence); ++i) {
            match = false;
            while (node != &TERM1 && node != &TERM0) {
                if (node->label == py::extract<char>(sequence[i])) {
                    node = node->branch1;
                    match = true;
                    break;
                }
                node = node->branch0;
            }
            if (!match) {
                return false;
            }
        }

        if (node != &TERM1) {
            return false;
        }

        return true;
    }
}; // namespace seqbdd

namespace seqbdd { namespace algorithm {
    typedef std::map<std::string, std::tuple<const Node*, int> > Failure;

    const Node* suffixdd(const py::str& sequence) {
        const Node* r = &TERM1;
        const Node* s = &TERM1;
        for (int i=py::len(sequence)-1; i >= 0; --i) {
            r = get_node(py::extract<char>(sequence[i]), r, &TERM1);
            s = union_(s, r);
        }

        return s;
    }

    const Node* next_node(const Node* node, const std::string& queue) {
        while (node != &TERM0 && node != &TERM1) {
            if (node->label == queue[0]) {
                if (queue.size() == 1) {
                    return node->branch1;
                } else {
                    return next_node(node->branch1, queue.substr(1));
                }
            }
            node = node->branch0;
        }

        return NULL;
    }

    void make(const Node* root, const Node* node, Failure& failure, std::string& queue) {
        bool found = false;
        for (size_t i=1; i < queue.size(); ++i) {
            const Node* next = next_node(root, queue.substr(i));
            if (next != NULL) {
                failure[queue] = std::tuple<const Node*, int>(next, i);
                found = true;
                break;
            }
        }
        if (!found) {
            failure[queue] = std::tuple<const Node*, int>(root, queue.size());
        }

        if (node != &TERM1) {
            if (node->branch0 != &TERM0) {
                make(root, node->branch0, failure, queue);
            }
            if (node->branch1 != &TERM0) {
                queue.push_back(node->label);
                make(root, node->branch1, failure, queue);
                queue.pop_back();
            }
        }
    }

    Failure make_failure(const Node* root) {
        Failure failure;
        std::string queue;

        make(root, root, failure, queue);

        return failure;
    }

    py::list search(const Node* root, Failure& failure, const std::string& sequence) {
        py::list result;

        const Node* node = root;
        std::string string = sequence + " ";
        std::string queue;
        size_t i = 0;
        while (i < string.size()) {
            // Check 0-branch to decide whether current node is terminal.
            const Node* node_b0 = node->branch0;
            while (node_b0 != &TERM0) {
                if (node_b0 == &TERM1) {
                    result.append(make_tuple(i-queue.size(), py::str(queue.c_str())));
                    break;
                }
                node_b0 = node_b0->branch0;
            }

            if (node->label == string[i]) {
                node = node->branch1;
                queue.push_back(string[i++]);
            } else {
                node = node->branch0;

                // Use failure function when label is not matched.
                if (node == &TERM0 || node == &TERM1) {
                    std::tuple<const Node*, int> tuple = failure[queue];
                    node = std::get<0>(tuple);
                    int n_pop = std::get<1>(tuple);

                    queue = queue.substr(n_pop);
                    i -= n_pop - 1;
                }
            }

            // Check terminal node.
            while (node == &TERM1) {
                result.append(make_tuple(i-queue.size(), py::str(queue.c_str())));

                std::tuple<const Node*, int> tuple = failure[queue];
                node = std::get<0>(tuple);
                int n_pop = std::get<1>(tuple);
                queue = queue.substr(n_pop);
                i -= n_pop - 1;
            }
        }

        return result;
    }
}}; // namespace alignment, seqbdd

namespace seqbdd { namespace alignment {
    // Cell value for direction tables.
    enum Cell {
        end = 0,
        match,
        insert,
        delete_,
    };

    using NodeVector = std::vector<const Node*>;
    using NodeIndex = std::map<const Node*, size_t>;
    using BestScore = std::map<std::string, int>;
    template<typename T> using Row =  std::vector<T>;
    template<typename T> using Matrix = std::vector<Row<T> >;

    // Gap symbol.
    static const char GAP = '-';

    template<typename T>
    inline bool has_not_item(std::vector<T> v, T x) {
        return std::find(v.begin(), v.end(), x) == v.end();
    }

    NodeVector get_nodes(const Node* root) {
        NodeVector result { &TERM0, &TERM1 };

        if (root == &TERM0 || root == &TERM1) {
            return result;
        }

        const Node* node = root;
        NodeVector stack { root };
        while (!stack.empty()) {
            node = stack.back();

            if (has_not_item(result, node->branch0)) {
                stack.push_back(node->branch0);
            } else if (has_not_item(result, node->branch1)) {
                stack.push_back(node->branch1);
            } else {
                result.push_back(node);
                stack.pop_back();
            }
        }

        return result;
    }

    inline NodeIndex get_node_index(const NodeVector& nodes) {
        NodeIndex node_index;

        for (size_t i=0; i < nodes.size(); ++i) {
            node_index[nodes[i]] = i;
        }

        return node_index;
    }

    inline int calc_score(const std::string& x, const std::string& y,
                          int gapopen, int gapext) {
        int score = 0;
        size_t length = std::min(x.size(), y.size());
        for (size_t i=0; i < length; ++i) {
            if (x[i] == GAP || y[i] == GAP) {
                if (i > 0 && (x[i-1] == GAP || y[i-1] == GAP)) {
                    score -= gapext;
                } else {
                    score -= gapopen;
                }
            } else {
                score += matrix::BLOSUM62::at(x[i], y[i]);
            }
        }

        return score;
    }

    BestScore get_best_score(const Node* root, int gapopen, int gapext) {
        BestScore best_score;
        std::vector<std::string> sequences = get_sequences(root);

        for (auto sequence : sequences) {
            best_score[sequence] = calc_score(sequence, sequence, gapopen, gapext);
        }

        return best_score;
    }

    template<typename T>
    inline Matrix<T> new_matrix(size_t height, size_t width) {
        return Matrix<T>(height, Row<T>(width));
    }

    inline void init_table(Matrix<int>& mat, size_t height, size_t width,
                           const NodeVector& nodes, NodeIndex& node_index,
                           int gapopen, int gapext) {
        for (size_t j=0; j < width; ++j) {
            mat[0][j] = 0;
        }

        for (size_t i=1; i < height; ++i) {
            size_t index = node_index[nodes[i+1]->branch1] - 1;
            if (index == 0) {
                mat[i][0] = -gapopen;
            } else {
                mat[i][0] = mat[index][0] - gapext;
            }
        }
    }

    inline void init_table_v(Matrix<int>& mat, size_t width, int gapopen, int gapext) {
        for (size_t j=0; j < width; ++j) {
            mat[0][j] = gapext - gapopen;
        }
    }

    inline void init_table_h(Matrix<int>& mat, const Matrix<int>& table,
                             size_t height, int gapopen, int gapext) {
        for (size_t i=1; i < height; ++i) {
            mat[i][0] = table[i][0] + gapext - gapopen;
        }
    }

    inline void init_dtable(Matrix<Cell>& mat, size_t height, size_t width) {
        for (size_t i=1; i < height; ++i) {
            mat[i][0] = Cell::delete_;
        }
        /*
        for (size_t j=0; j < width; ++j) {
            mat[0][j] = Cell::end;
        }
        */
    }

    inline void init_dtable_v(Matrix<Cell>& mat, size_t width) {
        /*
        for (size_t j=0; j < width; ++j) {
            mat[0][j] = Cell::end;
        }
        */
    }

    inline void init_dtable_h(Matrix<Cell>& mat, size_t height) {
        //mat[0][0] = Cell::end;
        for (size_t i=1; i < height; ++i) {
            mat[i][0] = Cell::delete_;
        }
    }

    inline void init_trace(Matrix<size_t>& mat, size_t height, size_t width,
                           const NodeVector& nodes, NodeIndex& node_index) {
        for (size_t i=1; i < height; ++i) {
            for (size_t j=0; j < width; ++j) {
                mat[i][j] = node_index[nodes[i+1]->branch1] - 1;
            }
        }
    }

    inline void fill_row(char vi, const std::string& sequence,
                         Row<int>& row, Row<int>& row_v, Row<int>& row_h,
                         Row<Cell>& drow, Row<Cell>& drow_v, Row<Cell>& drow_h,
                         Row<int>& b1_row, Row<int>& b1_row_v,
                         size_t width, int gapopen, int gapext) {
        for (size_t j=1; j < width; ++j) {
            int v_match = b1_row[j] - gapopen;
            int v_delete = b1_row_v[j] - gapext;
            if (v_match > v_delete) {
                row_v[j] = v_match;
                drow_v[j] = Cell::match;
            } else {
                row_v[j] = v_delete;
                drow_v[j] = Cell::delete_;
            }

            int h_match = row[j-1] - gapopen;
            int h_insert = row_h[j-1] - gapext;
            if (h_match > h_insert) {
                row_h[j] = h_match;
                drow_h[j] = Cell::match;
            } else {
                row_h[j] = h_insert;
                drow_h[j] = Cell::insert;
            }

            int match = b1_row[j-1] + matrix::BLOSUM62::at(vi, sequence[j-1]);
            int delete_ = row_v[j];
            int insert = row_h[j];
            if (match > insert) {
                if (match > delete_) {
                    row[j] = match;
                    drow[j] = Cell::match;
                } else {
                    row[j] = delete_;
                    drow[j] = Cell::delete_;
                }
            } else {
                if (insert > delete_) {
                    row[j] = insert;
                    drow[j] = Cell::insert;
                } else {
                    row[j] = delete_;
                    drow[j] = Cell::delete_;
                }
            }
        }
    }

    inline void fill_table(Matrix<int>& table, Matrix<int>& table_v, Matrix<int>& table_h,
                           Matrix<Cell>& dtable, Matrix<Cell>& dtable_v, Matrix<Cell>& dtable_h,
                           Matrix<size_t>& trace, size_t height, size_t width,
                           const NodeVector& nodes, NodeIndex& node_index,
                           const std::string& sequence, int gapopen, int gapext) {
        for (size_t i=1; i < height; ++i) {
            const Node* node = nodes[i+1];
            char vi = node->label;

            // 1-branch
            size_t b1_index = node_index[node->branch1] - 1;

            Row<int>& row = table[i];
            Row<int>& row_v = table_v[i];
            Row<int>& row_h = table_h[i];
            Row<Cell>& drow = dtable[i];
            Row<Cell>& drow_v = dtable_v[i];
            Row<Cell>& drow_h = dtable_h[i];

            fill_row(vi, sequence, row, row_v, row_h, drow, drow_v, drow_h,
                     table[b1_index], table_v[b1_index], width, gapopen, gapext);

            // 0-branch
            size_t b0_index = node_index[node->branch0];
            if (b0_index > 0) {
                b0_index--;

                Row<int>& b0_row = table[b0_index];
                Row<int>& b0_row_v = table_v[b0_index];
                Row<int>& b0_row_h = table_h[b0_index];
                Row<Cell>& b0_drow = dtable[b0_index];
                Row<Cell>& b0_drow_v = dtable_v[b0_index];
                Row<Cell>& b0_drow_h = dtable_h[b0_index];
                Row<size_t>& trace_row = trace[i];

                for (size_t j=0; j < width; ++j) {
                    if (b0_row[j] > row[j]) {
                        row[j] = b0_row[j];
                        row_v[j] = b0_row_v[j];
                        row_h[j] = b0_row_h[j];
                        drow[j] = b0_drow[j];
                        drow_v[j] = b0_drow_v[j];
                        drow_h[j] = b0_drow_h[j];

                        trace_row[j] = b0_index;
                    }
                }
            }
        }
    }

    inline void traceback(const Matrix<Cell>& dtable, const Matrix<Cell>& dtable_v,
                          const Matrix<Cell>& dtable_h, const Matrix<size_t>& trace,
                          size_t i, size_t j,
                          const NodeVector& nodes, NodeIndex& node_index,
                          const std::string& sequence,
                          std::string& seq_prime, std::string& query_prime) {
        const Matrix<Cell>* now = &dtable;

        query_prime.clear();
        seq_prime.clear();
        while ((*now)[i][j] != Cell::end) {
            // Jump row when replacing the row.
            while (i > 0 && trace[i][j] != node_index[nodes[i+1]->branch1]-1) {
                i = trace[i][j];
            }

            switch ((*now)[i][j]) {
            case Cell::match:
                if (now == &dtable) {
                    query_prime.push_back(nodes[i+1]->label);
                    seq_prime.push_back(sequence[j-1]);
                    i = trace[i][j];
                    --j;
                } else if (now == &dtable_v) {
                    query_prime.push_back(nodes[i+1]->label);
                    seq_prime.push_back(GAP);
                    i = trace[i][j];
                } else {
                    query_prime.push_back(GAP);
                    seq_prime.push_back(sequence[--j]);
                }

                now = &dtable;
                break;
            case Cell::delete_:
                if (dtable_v[i][j] == Cell::delete_) {
                    now = &dtable_v;
                }

                query_prime.push_back(nodes[i+1]->label);
                seq_prime.push_back(GAP);
                i = trace[i][j];
                break;
            case Cell::insert:
                if (dtable_h[i][j] == Cell::insert) {
                    now = &dtable_h;
                }

                query_prime.push_back(GAP);
                seq_prime.push_back(sequence[--j]);
                break;
            case Cell::end:
                break;
            }
        }
    }

    inline int get_floor_threshold(const BestScore& best_score, float percent) {
        auto min = *std::min_element(best_score.begin(), best_score.end(),
                                     [](const std::pair<std::string, int>& x,
                                        const std::pair<std::string, int>& y) {
                                        return x.second < y.second; });
        return static_cast<int>(min.second * percent);
    }

    inline std::vector<size_t> start_points(const Row<int>& row, size_t width) {
        std::vector<size_t> starts(width);

        for (size_t j=0; j < width; ++j) {
            starts[j] = j;
        }
        std::sort(starts.begin(), starts.end(),
                  [&](size_t x, size_t y) -> bool { return row[x] > row[y]; });

        return starts;
    }

    py::list glocal_alignment(const std::string& sequence, const Node* query,
                              const NodeVector& nodes, BestScore& best_score,
                              int gapopen, int gapext, float p) {
        NodeIndex node_index = get_node_index(nodes);

        size_t height = nodes.size() - 1;
        size_t width = sequence.size() + 1;

        // Score tables
        Matrix<int> table = new_matrix<int>(height, width);
        Matrix<int> table_v = new_matrix<int>(height, width);
        Matrix<int> table_h = new_matrix<int>(height, width);
        init_table(table, height, width, nodes, node_index, gapopen, gapext);
        init_table_v(table_v, width, gapopen, gapext);
        init_table_h(table_v, table, height, gapopen, gapext);

        // Direction tables
        Matrix<Cell> dtable = new_matrix<Cell>(height, width);
        Matrix<Cell> dtable_v = new_matrix<Cell>(height, width);
        Matrix<Cell> dtable_h = new_matrix<Cell>(height, width);
        init_dtable(dtable, height, width);
        init_dtable_v(dtable, width);
        init_dtable_h(dtable, height);

        // Trace table
        Matrix<size_t> trace = new_matrix<size_t>(height, width);
        init_trace(trace, height, width, nodes, node_index);

        std::string rev_seq(sequence);
        std::reverse(rev_seq.begin(), rev_seq.end());
        fill_table(table, table_v, table_h, dtable, dtable_v, dtable_h,
                   trace, height, width, nodes, node_index,
                   rev_seq, gapopen, gapext);

        py::list results;
        std::string seq_prime;
        std::string query_prime;

        int floor = get_floor_threshold(best_score, p);

        int i = height - 1;
        std::vector<size_t> starts = start_points(table[i], width);
        for (size_t j : starts) {
            int score = table[i][j];
            if (score < floor) {
                break;
            }

            traceback(dtable, dtable_v, dtable_h, trace,
                      i, j, nodes, node_index, rev_seq,
                      seq_prime, query_prime);

            // Remove gaps
            std::string raw_query;
            std::string raw_seq;
            std::copy_if(query_prime.begin(), query_prime.end(),
                         std::back_inserter(raw_query),
                         [](char x) { return x != GAP; });
            std::copy_if(seq_prime.begin(), seq_prime.end(),
                         std::back_inserter(raw_seq),
                         [](char x) { return x != GAP; });

            if (score >= best_score[raw_query]*p) {
                int seq_index = sequence.size() - j;

                results.append(
                    py::make_tuple(raw_query,
                                   sequence,
                                   query_prime,
                                   seq_prime,
                                   py::make_tuple(0, raw_query.size()),
                                   py::make_tuple(seq_index, seq_index+raw_seq.size()),
                                   score));
            }
        }

        return results;
    }
}}; //namespace alignment, seqbdd

BOOST_PYTHON_MODULE(_seqbdd) {
    using namespace seqbdd;

    // Types
    py::class_<Node>("Node", py::init<char, const Node*, const Node*>());

    py::class_<algorithm::Failure>("Failure")
        .def(py::map_indexing_suite<algorithm::Failure>());
    py::class_<alignment::BestScore>("BestScore")
        .def(py::map_indexing_suite<alignment::BestScore>());
    py::class_<alignment::NodeVector>("NodeVector")
        .def(py::vector_indexing_suite<alignment::NodeVector>());

    // Interface by py::list
    py::def("from_sequences", from_sequences,
            py::return_value_policy<py::reference_existing_object>());
    py::def("values", values,
            py::return_value_policy<py::return_by_value>());

    // Set operations
    py::def("union_", union_,
            py::return_value_policy<py::reference_existing_object>());
    py::def("intersection", intersection,
            py::return_value_policy<py::reference_existing_object>());
    py::def("difference", difference,
            py::return_value_policy<py::reference_existing_object>());
    py::def("symmetric_difference", symmetric_difference,
            py::return_value_policy<py::reference_existing_object>());
    py::def("one_terminal", one_terminal,
            py::return_value_policy<py::reference_existing_object>());
    py::def("zero_terminal", zero_terminal,
            py::return_value_policy<py::reference_existing_object>());
    py::def("onset", onset,
            py::return_value_policy<py::reference_existing_object>());
    py::def("push", push,
            py::return_value_policy<py::reference_existing_object>());
    py::def("top", top,
            py::return_value_policy<py::return_by_value>());
    py::def("count", count,
            py::return_value_policy<py::return_by_value>());
    py::def("count_node", count_node,
            py::return_value_policy<py::return_by_value>());

    // Implements Python statements
    py::def("has_sequence", has_sequence,
            py::return_value_policy<py::return_by_value>());
    py::def("equal_nodes", equal_nodes,
            py::return_value_policy<py::return_by_value>());
    py::def("not_equal_nodes", not_equal_nodes,
            py::return_value_policy<py::return_by_value>());

    // Algorithms
    // [SuffixDD]
    py::def("suffixdd", algorithm::suffixdd,
            py::return_value_policy<py::reference_existing_object>());

    // [Aho-Corasick for SeqBDDs]
    py::def("make_failure", algorithm::make_failure,
            py::return_value_policy<py::return_by_value>());
    py::def("search", algorithm::search,
            py::return_value_policy<py::return_by_value>());

    // [Glocal sequence alignment using SeqBDDs]
    py::def("get_nodes", alignment::get_nodes,
            py::return_value_policy<py::return_by_value>());
    py::def("get_best_score", alignment::get_best_score,
            py::return_value_policy<py::return_by_value>());
    py::def("glocal_alignment", alignment::glocal_alignment,
            py::return_value_policy<py::return_by_value>());
}

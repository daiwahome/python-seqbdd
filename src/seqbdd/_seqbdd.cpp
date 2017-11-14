#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

namespace py = boost::python;

typedef struct _Node {
    char label;
    const struct _Node* branch1;
    const struct _Node* branch0;

    _Node(char l, const struct _Node* b1, const struct _Node* b0)
        : label(l), branch1(b1), branch0(b0) {}
} Node;

static const Node TERM0('0', NULL, NULL);
static const Node TERM1('1', NULL, NULL);

class nodehash {
public:
    size_t operator()(const std::tuple<char, const Node*, const Node*>& x) const {
        return (size_t)std::get<1>(x) ^ (size_t)std::get<2>(x);
    }
};

static std::unordered_map<std::tuple<char, const Node*, const Node*>,
    std::unique_ptr<const Node>, nodehash> cache;

const Node* get_node(char x, const Node* p1, const Node* p0) {
    if (p1 == &TERM0) {
        return p0;
    } else {
        std::tuple<char, const Node*, const Node*> key(x, p1, p0);

        if (!cache.count(key)) {
            cache[key] = std::unique_ptr<const Node>(new Node(x, p1, p0));
        }

        return cache[key].get();
    }
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

void trace(const Node* node, py::list& sequences, std::vector<char>& queue) {
    if (node == &TERM1) {
        sequences.append(py::str(std::string(queue.begin(), queue.end()).c_str()));
    } else if (node != &TERM0) {
        trace(node->branch0, sequences, queue);

        queue.push_back(node->label);
        trace(node->branch1, sequences, queue);
        queue.pop_back();
    }
}

py::list values(const Node* root) {
    py::list sequences;
    std::vector<char> queue;

    trace(root, sequences, queue);

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

void make(const Node* root, const Node* node,
          std::map<std::string, std::tuple<const Node*, int> >& failure,
          std::string& queue) {
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

std::map<std::string, std::tuple<const Node*, int> > make_failure(const Node* root) {
    std::map<std::string, std::tuple<const Node*, int> > failure;
    std::string queue;

    make(root, root, failure, queue);

    return failure;
}

py::list search(const Node* root,
                std::map<std::string, std::tuple<const Node*, int> >& failure,
                const std::string& sequence) {
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

BOOST_PYTHON_MODULE(_seqbdd) {
    py::class_<Node>("Node", py::init<char, const Node*, const Node*>());
    py::class_<std::map<std::string, std::tuple<const Node*, int> > >("failure")
        .def(py::map_indexing_suite<std::map<std::string, std::tuple<const Node*, int> > >());

    py::def("union_", union_,
            py::return_value_policy<py::reference_existing_object>());
    py::def("intersection", intersection,
            py::return_value_policy<py::reference_existing_object>());
    py::def("difference", difference,
            py::return_value_policy<py::reference_existing_object>());
    py::def("symmetric_difference", symmetric_difference,
            py::return_value_policy<py::reference_existing_object>());
    py::def("from_sequences", from_sequences,
            py::return_value_policy<py::reference_existing_object>());
    py::def("one_terminal", one_terminal,
            py::return_value_policy<py::reference_existing_object>());
    py::def("zero_terminal", zero_terminal,
            py::return_value_policy<py::reference_existing_object>());
    py::def("onset", onset,
            py::return_value_policy<py::reference_existing_object>());
    py::def("push", push,
            py::return_value_policy<py::reference_existing_object>());
    py::def("suffixdd", suffixdd,
            py::return_value_policy<py::reference_existing_object>());
    py::def("make_failure", make_failure,
            py::return_value_policy<py::return_by_value>());
    py::def("search", search,
            py::return_value_policy<py::return_by_value>());
    py::def("top", top,
            py::return_value_policy<py::return_by_value>());
    py::def("count", count,
            py::return_value_policy<py::return_by_value>());
    py::def("count_node", count_node,
            py::return_value_policy<py::return_by_value>());
    py::def("values", values,
            py::return_value_policy<py::return_by_value>());
    py::def("has_sequence", has_sequence,
            py::return_value_policy<py::return_by_value>());
    py::def("equal_nodes", equal_nodes,
            py::return_value_policy<py::return_by_value>());
    py::def("not_equal_nodes", not_equal_nodes,
            py::return_value_policy<py::return_by_value>());
}

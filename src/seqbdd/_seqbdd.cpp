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
    struct _Node* branch1;
    struct _Node* branch0;

    _Node(char l, struct _Node* b1, struct _Node* b0) {
        label = l;
        branch1 = b1;
        branch0 = b0;
    }
} Node;

static Node TERM0('0', NULL, NULL);
static Node TERM1('1', NULL, NULL);

class nodehash {
public:
    size_t operator()(const std::tuple<char, Node*, Node*>& x) const {
        return (size_t)std::get<1>(x) ^ (size_t)std::get<2>(x);
    }
};

static std::unordered_map<std::tuple<char, Node*, Node*>, std::shared_ptr<Node>, nodehash> cache;

Node* get_node(char x, Node* p1, Node* p0) {
    if (p1 == &TERM0) {
        return p0;
    } else {
        std::tuple<char, Node*, Node*> key(x, p1, p0);

        if (!cache.count(key)) {
            cache[key] = std::shared_ptr<Node>(new Node(x, p1, p0));
        }

        return cache[key].get();
    }
}

int sign(Node* x) {
    if (x == &TERM0) {
        return 0;
    } else {
        return 1;
    }
}

template<typename Func>
Node* meld(Func op, Node* p, Node* q) {
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

Node* union_(Node* p, Node* q) {
    return meld([](int x, int y) -> int { return x+y; }, p, q);
}

Node* intersection(Node* p, Node* q) {
    return meld([](int x, int y) -> int { return x*y; }, p, q);
}

Node* difference(Node* p, Node* q) {
    return meld([](int x, int y) -> int { return x&(!y); }, p, q);
}

Node* symmetric_difference(Node* p, Node* q) {
    return meld([](int x, int y) -> int { return x^y; }, p, q);
}

Node* one_terminal() {
    return &TERM1;
}

Node* zero_terminal() {
    return &TERM0;
}

Node* onset(Node* root, char x) {
    Node* node = root;
    while (node != &TERM1 && node != &TERM0) {
        if (node->label == x) {
            return node->branch1;
        }
        node = node->branch0;
    }

    return node;
}

Node* push(Node* root, char x) {
    return get_node(x, root, &TERM0);
}

char top(Node* root) {
    return root->label;
}

int count(Node* root) {
    if (root == &TERM0) {
        return 0;
    } else if (root == &TERM1) {
        return 1;
    }

    Node* node = root;
    std::unordered_map<Node*, int> counter { { &TERM0, 0 }, { &TERM1, 1 } };
    std::vector<Node*> stack { root };
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

int count_node(Node* root) {
    if (root == &TERM0 || root == &TERM1) {
        return 0;
    }

    Node* node = root;
    std::unordered_set<Node*> set { &TERM0, &TERM1 };
    std::vector<Node*> stack { root };
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

Node* make_string(const py::str& sequence) {
    Node* root = &TERM1;
    for (int i=py::len(sequence)-1; i >= 0; --i) {
        root = get_node(py::extract<char>(sequence[i]), root, &TERM0);
    }

    return root;
}

Node* from_sequences(const py::list& sequences) {
    Node* root = &TERM0;
    for (int i=0; i < py::len(sequences); ++i) {
        root = union_(root, make_string(py::extract<py::str>(sequences[i])));
    }

    return root;
}

void trace(Node* node, py::list& sequences, std::vector<char>& queue) {
    if (node->label == '1') {
        sequences.append(py::str(std::string(queue.begin(), queue.end()).c_str()));
    } else if (node->label != '0') {
        trace(node->branch0, sequences, queue);

        queue.push_back(node->label);
        trace(node->branch1, sequences, queue);
        queue.pop_back();
    }
}

py::list values(Node* root) {
    py::list sequences;
    std::vector<char> queue;

    trace(root, sequences, queue);

    return sequences;
}

bool has_sequence(Node* root, const py::str& sequence) {
    bool match = false;

    Node* node = root;
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

Node* suffixdd(const py::str& sequence) {
    Node* r = &TERM1;
    Node* s = &TERM1;
    for (int i=py::len(sequence)-1; i >= 0; --i) {
        r = get_node(py::extract<char>(sequence[i]), r, &TERM1);
        s = union_(s, r);
    }

    return s;
}

Node* next_node(Node* node, const std::string& queue) {
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

void make(Node* root, Node* node,
          std::map<std::string, std::tuple<Node*, int> >& failure,
          std::string& queue) {
    bool found = false;
    for (size_t i=1; i < queue.size(); ++i) {
        Node* next = next_node(root, queue.substr(i));
        if (next != NULL) {
            failure[queue] = std::tuple<Node*, int>(next, i);
            found = true;
            break;
        }
    }
    if (!found) {
        failure[queue] = std::tuple<Node*, int>(root, queue.size());
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

std::map<std::string, std::tuple<Node*, int> > make_failure(Node* root) {
    std::map<std::string, std::tuple<Node*, int> > failure;
    std::string queue;

    make(root, root, failure, queue);

    return failure;
}

py::list search(Node* root,
                std::map<std::string, std::tuple<Node*, int> >& failure,
                const std::string& sequence) {
    py::list result;

    Node* node = root;
    std::string string = sequence + " ";
    std::string queue;
    size_t i = 0;
    while (i < string.size()) {
        // Check 0-branch to decide whether current node is terminal.
        Node* node_b0 = node->branch0;
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
                std::tuple<Node*, int> tuple = failure[queue];
                node = std::get<0>(tuple);
                int n_pop = std::get<1>(tuple);

                queue = queue.substr(n_pop);
                i -= n_pop - 1;
            }
        }

        // Check terminal node.
        while (node == &TERM1) {
            result.append(make_tuple(i-queue.size(), py::str(queue.c_str())));

            std::tuple<Node*, int> tuple = failure[queue];
            node = std::get<0>(tuple);
            int n_pop = std::get<1>(tuple);
            queue = queue.substr(n_pop);
            i -= n_pop - 1;
        }
    }

    return result;
}

BOOST_PYTHON_MODULE(_seqbdd) {
    py::class_<Node>("Node", py::init<char, Node*, Node*>());
    py::class_<std::map<std::string, std::tuple<Node*, int> > >("failure")
        .def(py::map_indexing_suite<std::map<std::string, std::tuple<Node*, int> > >());

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
}
